import logging

import pandas as pd
import pyfocus as pf

from sqlalchemy import create_engine


__all__ = ["import_fusion", "import_predixcan"]

# TODO: implement exporting to predixcan/fusion


def import_fusion(path, name, tissue, assay, session, bulk_amt=100):
    """
    Import weights from a PrediXcan db into the FOCUS db.

    :param path:  string path to the PrediXcan db
    :param tissue: str name of the tissue
    :param assay: technology assay to measure abundance
    :param session: sqlalchemy.Session object for the FOCUS db
    :param bulk_amt: int the number of objects to store until bulk commit

    :return:  None
    """
    import re
    import os
    import warnings

    try:
        import mygene
        # suppress warnings about R build
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import rpy2.robjects as robj
    except ImportError:
        log.error("Import submodule requires mygene and rpy2 to be installed.")
        raise

    import numpy as np


    log = logging.getLogger(pf.LOG)

    log.info("Starting import from FUSION database {}".format(path))
    db_ref_panel = pf.RefPanel(ref_name=name, tissue=tissue, assay=assay)
    ses = None

    load_func = robj.r['load']
    local_dir = os.path.dirname(os.path.abspath(path))

    # we need this to grab Ensembl IDs for genes
    mg = mygene.MyGeneInfo()

    db_objs = []

    with open(path, "rt") as fusion_info:
        header = fusion_info.readline()

        for line in fusion_info:
            wgt_name, g_name, chrom, txstart, txstop = line.split()

            # WGT ID CHR P0 P1
            # METSIM.ADIPOSE.RNASEQ/METSIM.LINC00115.wgt.RDat LINC00115 1 761586 762902
            log.debug("Importing {} model".format(wgt_name))

            # this call should create the following:
            # 'wgt.matrix', 'snps', 'cv.performance', 'hsq', and 'hsq.pv'
            wgt_path = "{}/{}".format(local_dir, wgt_name)

            # this is more of a bottleneck than querying the web using mygene
            # if it is more IO than computational we can push this to prod/consum queue + thread approach
            # but likely not worth it since importing is done rarely
            load_func(wgt_path)

            # this mygene API is so nice...
            result = mg.query(g_name, scopes='symbol', fields=['ensembl.gene,genomic_pos,symbol,ensembl.type_of_gene,alias'],
                              species="human")

            gene_info = dict()
            id_dict = dict()
            # hits are ordered by match quality.
            for hit in result['hits']:
                if "ensembl" not in hit:
                    # nothing in db
                    continue

                if hit["symbol"] != g_name and "alias" in hit and g_name not in hit["alias"]:
                    # not direct match
                    continue

                if "genomic_pos" not in hit:
                    continue

                ens = hit["ensembl"]
                pos = hit["genomic_pos"]

                # sometimes we have multiple ENSG entries due to diff haplotypes.
                # just reduce the single-case to the multi by a singleton list
                if type(ens) is dict:
                    ens = [ens]
                if type(pos) is dict:
                    pos = [pos]

                # grab the type for when we match against pos
                for e_hit in ens:
                    id_dict[e_hit["gene"]] = e_hit["type_of_gene"]

                for p_hit in pos:
                    if not re.match("[0-9]{1,2}|X|Y", p_hit["chr"], re.IGNORECASE):
                        continue
                    g_id = p_hit["ensemblgene"]
                    g_type = id_dict.get(p_hit["ensemblgene"])

                    if len(gene_info) == 0:
                        # grab any info if we haven't yet
                        gene_info["geneid"] = g_id
                        gene_info["type"] = g_type
                    elif "protein" in g_type:
                        # prioritize protein coding and break out if we find one
                        gene_info["geneid"] = g_id
                        gene_info["type"] = g_type

            if len(gene_info) == 0:
                # we didn't get any hits from our query
                # just use the gene-name as ens-id...
                log.warning("Unable to match {} to Ensembl ID. Using symbol for ID".format(g_name))
                gene_info["geneid"] = g_name
                gene_info["type"] = None

            # build info on the gene
            gene_info["txid"] = None
            gene_info["name"] = g_name
            gene_info["chrom"] = chrom
            gene_info["txstart"] = txstart
            gene_info["txstop"] = txstop

            # get the multi-SNP method with the best cvR2
            methods = list(robj.r['cv.performance'].colnames)
            types = list(robj.r['cv.performance'].rownames)
            if "rsq" not in types:
                raise ValueError("No R2 value for model {}".format(path))
            if "pval" not in types:
                raise ValueError("No R2 p-value for model {}".format(path))

            rsq_idx = types.index("rsq")
            pval_idx = types.index("pval")
            values = np.array(robj.r['cv.performance'])

            r2 = 0
            method = None
            r2idx = 0
            for idx, value in enumerate(values[rsq_idx]):
                if methods[idx] == "top1":
                    continue
                if value > r2:
                    r2 = value
                    method = methods[idx]
                    r2idx = idx
                    pval = values[pval_idx, idx]

            # keep attributes
            attrs = dict()
            attrs["cv.R2"] = r2
            attrs["cv.R2.pval"] = pval

            # grab the actual weights
            wgts = np.array(robj.r['wgt.matrix']).T[r2idx]

            # SNPs data frame
            # V1 V2 V3 V4 V5 V6
            # 11 rs2729762 0 77033699 G A
            snps = robj.r['snps']
            snp_info = pd.DataFrame({"snp": list(snps[1]),
                                     "chrom": list(snps[0]),
                                     "pos": list(snps[3]),
                                     "a1": list(snps[4]),
                                     "a0": list(snps[5])})

            # if we're using a sparse model there is no need to store info on zero'd SNPs
            keep = np.logical_not(np.isclose(wgts, 0))
            wgts = wgts[keep]
            snp_info = snp_info[keep]

            model = pf.build_model(gene_info, snp_info, db_ref_panel, wgts, ses, attrs, method)
            db_objs.append(model)
            if len(db_objs) == bulk_amt:
                session.bulk_save_objects(db_objs)
                try:
                    session.commit()
                    db_objs = []
                except Exception as comm_err:
                    session.rollback()
                    raise Exception("Failed committing to db")

    # get any remaining
    try:
        session.bulk_save_objects(db_objs)
        session.commit()
    except Exception as comm_err:
        session.rollback()

    log.info("Finished import from FUSION database {}".format(path))
    return


def export_fusion(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_fusion not implemented!")
    return


def import_predixcan(path, name, tissue, assay, session, bulk_amt=100):
    """
    Import weights from a PrediXcan db into the FOCUS db.

    :param path:  string path to the PrediXcan db
    :param tissue: str name of the tissue
    :param assay: technology assay to measure abundance
    :param session: sqlalchemy.Session object for the FOCUS db
    :param bulk_amt: int the number of objects to store until bulk commit

    :return:  None
    """
    log = logging.getLogger(pf.LOG)

    import re
    import numpy as np
    try:
        import mygene
    except ImportError:
        log.error("Import submodule requires mygene and rpy2 to be installed.")
        raise

    log.info("Starting import from PrediXcan database {}".format(path))
    pred_engine = create_engine("sqlite:///{}".format(path))

    weights = pd.read_sql_table('weights', pred_engine)
    extra = pd.read_sql_table('extra', pred_engine)

    mg = mygene.MyGeneInfo()

    db_objs = []

    db_ref_panel = pf.RefPanel(ref_name=name, tissue=tissue, assay=assay)
    ses = None
    method = "ElasticNet"

    for gid, gene in weights.groupby("gene"):
        log.debug("Importing gene model {}".format(gid))
        gene_extra = extra.loc[extra.gene == gid]

        chrom = gene.varID.values[0].split("_")[0]  # grab chromosome from varID
        pos = gene.varID.map(lambda x: int(x.split("_")[1])).values  # grab basepair pos

        g_id = gene_extra.gene.values[0]
        g_name = gene_extra.genename.values[0]
        p_idx = g_id.rfind(".")
        if p_idx > -1:
            query_id = g_id[:p_idx]
        else:
            query_id = g_id

        txstart = txstop = np.median(pos)
        result = mg.query(query_id, scopes="ensembl.gene", fields=["genomic_pos_hg19,symbol,alias"], species="human")
        for hit in result['hits']:
            if hit["symbol"] != g_name and "alias" in hit and g_name not in hit["alias"]:
                continue
            if "genomic_pos_hg19" not in hit:
                continue

            gpos = hit["genomic_pos_hg19"]
            if type(gpos) is dict:
                gpos = [gpos]

            for entry in gpos:
                # skip non-primary assembles. they have weird CHR entries e.g., CHR_HSCHR1_1_CTG3
                if not re.match("[0-9]{1,2}|X|Y", entry["chr"], re.IGNORECASE):
                    continue

                txstart = entry['start']
                txstop = entry['end']
                break

            if txstart is not None:
                # we want to use standardized Ensembl identifiers; not GENCODE modified ones...
                g_id = query_id
                break

        gene_info = dict()
        gene_info["geneid"] = g_id
        gene_info["txid"] = None
        gene_info["name"] = g_name
        gene_info["type"] = gene_extra.gene_type.values[0]
        gene_info["chrom"] = chrom
        gene_info["txstart"] = txstart
        gene_info["txstop"] = txstop

        snp_info = pd.DataFrame({"snp": gene.rsid.values,
                                "chrom": [chrom] * len(gene),
                                "pos": pos,
                                "a1": gene.eff_allele.values,
                                "a0": gene.ref_allele.values})

        wgts = gene.weight.values

        attrs = dict()
        attrs["cv.R2"] = gene_extra["cv_R2_avg"].values[0]
        attrs["cv.R2.pval"] = gene_extra["nested_cv_fisher_pval"].values[0]

        # build model
        model = pf.build_model(gene_info, snp_info, db_ref_panel, wgts, ses, attrs, method)
        db_objs.append(model)
        if len(db_objs) == bulk_amt:
            session.bulk_save_objects(db_objs)
            try:
                session.commit()
                db_objs = []
            except Exception as comm_err:
                session.rollback()
                raise Exception("Failed committing to db")

    # get any remaining
    try:
        session.bulk_save_objects(db_objs)
        session.commit()
    except Exception as comm_err:
        session.rollback()

    log.info("Finished import from PrediXcan database {}".format(path))
    return


def export_predixcan(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_predixcan not implemented!")
    return
