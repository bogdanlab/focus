__all__ = ["is_file", "write_output"]


def is_file(f):
    """
    Check if f is a file-like object. Compatible with Python 2 and 3

    See https://stackoverflow.com/questions/47518707/check-if-an-object-is-a-file-in-python-2-and-3

    :param f:
    :return: True if `f` is file-like; False otherwise.
    """
    import sys
    return isinstance(f, file) if sys.version_info[0] == 2 else hasattr(f, 'read')


def write_output(results, output, append=False):
    results.to_csv(output, sep="\t", mode="a" if append else "w", header=not append, index=False)
    return
