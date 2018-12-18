import logging

import pyfocus as pf

__all__ = ["import_fusion", "import_metaxcan"]


def import_fusion(path, session):
    log = logging.getLogger(pf.LOG)
    try:
        import rpy2
    except ImportError as ie:
        raise ImportError("Import FUSION weights requires rpy2 library")

    from rpy2.robjects import pandas2ri as p2r
    r = rpy2.robjects.r

    load = r['load']
    with open(path, "r") as fusion_db:
        for line in fusion_db:
            load(line)

            snps = p2r.ri2py(r['snps'])


    pass

def import_metaxcan(path, session):
    log = logging.getLogger(pf.LOG)
    pass
