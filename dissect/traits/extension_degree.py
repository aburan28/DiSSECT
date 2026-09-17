from typing import List

from dissect.traits import Trait


def factor_degree(degree):
    """Prime factorization of a small positive integer, ascending with multiplicity."""
    factors = []
    remaining, divisor = degree, 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            factors.append(divisor)
            remaining //= divisor
        divisor += 1
    if remaining > 1:
        factors.append(remaining)
    return factors


def proper_divisors(degree):
    """Divisors d of `degree` with 1 < d < degree, ascending."""
    return [d for d in range(2, degree) if degree % d == 0]


class ExtensionDegreeTrait(Trait):
    NAME = "extension_degree"
    DESCRIPTION = (
        "Factorization of the degree $m$ of the field over its prime field. A composite $m$"
        " admits proper subfields, which is the precondition for a Weil descent (GHS) transfer"
        " of the discrete logarithm to a hyperelliptic curve over a subfield."
    )
    INPUT = {}
    OUTPUT = {
        "degree": (int, "Degree $m$ of the field over its prime field"),
        "factorization": (List[int], "Prime factorization of $m$"),
        "proper_divisor_count": (int, "Number of divisors $d$ with $1 < d < m$"),
        "largest_proper_divisor": (
            int,
            "Largest such divisor, or 0 when $m$ is 1 or prime",
        ),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Degree of the field extension and its divisor structure.

        A prime field has degree 1 and a prime-degree binary field has no proper
        divisor, so both report 0 for `largest_proper_divisor`. Any nonzero value
        names the subfield a descent would target.
        """
        degree = int(curve.field().degree())
        divisors = proper_divisors(degree)
        return {
            "degree": degree,
            "factorization": factor_degree(degree),
            "proper_divisor_count": len(divisors),
            "largest_proper_divisor": divisors[-1] if divisors else 0,
        }


def test_factor_degree():
    assert factor_degree(1) == []
    assert factor_degree(163) == [163]
    assert factor_degree(176) == [2, 2, 2, 2, 11]
    assert factor_degree(155) == [5, 31]


def test_proper_divisors():
    assert proper_divisors(1) == []
    assert proper_divisors(163) == []
    assert proper_divisors(155) == [5, 31]
    assert proper_divisors(176) == [2, 4, 8, 11, 16, 22, 44, 88]


def test_prime_degree_curves_have_no_descent_target():
    for degree in (113, 163, 233, 283, 409, 571):
        assert proper_divisors(degree) == []


def test_composite_degree_curves_do():
    # The X9.62 c2pnb family and the Oakley groups.
    for degree in (155, 176, 185, 208, 272, 304, 368):
        assert proper_divisors(degree) != []
