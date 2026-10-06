"""No-op numba stand-in for machines without numba (e.g. the D:\ machine, Python 3.11).

Use only via PYTHONPATH:  export PYTHONPATH="D:/Chem-Research/Np/tools/numba_stub"
@njit / @njit(...) return the plain Python function: identical results, slower. The project's
relativity.py / ks_atom.py numerics are cached (models/first_principles/cache), so normal use is fast.
"""


def njit(*args, **kwargs):
    if len(args) == 1 and callable(args[0]) and not kwargs:
        return args[0]
    return lambda f: f


jit = njit
prange = range
