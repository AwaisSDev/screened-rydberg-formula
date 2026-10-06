"""Angular coefficients for the Coulomb interaction (Condon-Shortley conventions).

    c^k(l m, l' m') = (-1)^m sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0) (l k l'; -m m-m' m')

The two-electron Coulomb matrix element between spin-orbitals (n l m sigma) is

    <ab|1/r12|cd> = delta(s_a,s_c) delta(s_b,s_d) delta(m_a+m_b, m_c+m_d)
                    sum_k c^k(l_a m_a, l_c m_c) c^k(l_d m_d, l_b m_b) R^k(ab; cd).

c^k(lm,lm) and c^k(lm,l'm')^2 are rational numbers; we return them exactly (Fraction) for the
diagonal (determinant) energies, and as floats for configuration-interaction matrices.
"""
from fractions import Fraction
from functools import lru_cache

import sympy
from sympy.physics.wigner import wigner_3j


@lru_cache(maxsize=None)
def ck_sym(l, m, l2, m2, k):
    if (l + k + l2) % 2 or k < abs(l - l2) or k > l + l2 or abs(m - m2) > k:
        return sympy.Integer(0)
    v = (-1) ** m * sympy.sqrt((2 * l + 1) * (2 * l2 + 1)) * wigner_3j(l, k, l2, 0, 0, 0) \
        * wigner_3j(l, k, l2, -m, m - m2, m2)
    return sympy.nsimplify(sympy.simplify(v))


@lru_cache(maxsize=None)
def ck_float(l, m, l2, m2, k):
    return float(ck_sym(l, m, l2, m2, k))


def _to_frac(x):
    x = sympy.nsimplify(sympy.simplify(x))
    assert x.is_Rational, x
    return Fraction(int(x.p), int(x.q))


@lru_cache(maxsize=None)
def ck_diag(l, m, k):
    """c^k(lm, lm), exact."""
    return _to_frac(ck_sym(l, m, l, m, k))


@lru_cache(maxsize=None)
def ck_sq(l, m, l2, m2, k):
    """c^k(lm, l'm')^2, exact."""
    return _to_frac(ck_sym(l, m, l2, m2, k) ** 2)


@lru_cache(maxsize=None)
def threej0_sq(l, k, l2):
    """(l k l'; 0 0 0)^2, exact."""
    return _to_frac(wigner_3j(l, k, l2, 0, 0, 0) ** 2)
