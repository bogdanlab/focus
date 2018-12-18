__all__ = ["write_output"]


def write_output(imputed_gwas, output, append=False):
    imputed_gwas.to_csv(output, sep="\t", mode="a" if append else "w", header=not append, index=False)
    return
