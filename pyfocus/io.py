__all__ = ["is_file", "write_output"]


def is_file(f):
    """
    Check if f is a file-like object. Compatible with Python 3

    :param f:
    :return: True if `f` is file-like; False otherwise.
    """
    return hasattr(f, 'read')


def write_output(results, output, append=False):
    results.to_csv(output, sep="\t", mode="a" if append else "w", header=not append, index=False, na_rep="NA",
                   float_format="%.3g")
    return
