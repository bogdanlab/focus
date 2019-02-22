import logging

import pandas as pd
import pyfocus as pf

from sqlalchemy import create_engine


__all__ = ["import_predixcan"]

# TODO: implement exporting to predixcan/fusion


def import_fusion(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("import_fusion not implemented!")
    return


def export_fusion(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_fusion not implemented!")
    return


def import_predixcan(path, name, tissue, assay, session):
    """
    Import weights from a PrediXcan db into the FOCUS db.

    :param path:  string path to the PrediXcan db
    :param tissue: str name of the tissue
    :param assay: technology assay to measure abundance
    :param session: sqlalchemy.Session object for the FOCUS db

    :return:  None
    """
    log = logging.getLogger(pf.LOG)

    log.info("Starting import from PrediXcan database {}".format(path))
    pred_engine = create_engine("sqlite:///{}".format(path))

    weights = pd.read_sql_table('weights', pred_engine)
    extra = pd.read_sql_table('extra', pred_engine)

    db_ref_panel = pf.RefPanel(ref_name=name, tissue=tissue, assay=assay)
    ses = None
    method = "ElasticNet"

    for gid, gene in weights.groupby("gene"):
        log.debug("Importing gene model {}".format(gid))
        gene_extra = extra.loc[extra.gene == gid]

        chrom = gene.varID.values[0].split("_")[0]  # grab chromosome from varID
        pos = gene.varID.map(lambda x: int(x.split("_")[1])).values  # grab basepair pos

        gene_info = dict()
        gene_info["geneid"] = gene_extra.gene.values[0]
        gene_info["txid"] = None
        gene_info["name"] = gene_extra.genename.values[0]
        gene_info["type"] = gene_extra.gene_type.values[0]
        gene_info["chrom"] = chrom
        gene_info["txstart"] = None
        gene_info["txstop"] = None

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
        session.add(model)
        try:
            session.commit()
        except Exception as comm_err:
            session.rollback()
            raise Exception("Failed committing to db")

    log.info("Finished import from PrediXcan database {}".format(path))
    return


def export_predixcan(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_predixcan not implemented!")
    return
