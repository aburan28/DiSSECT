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
    """Degrees of all proper subfields, including the prime field (d = 1)."""
    return [d for d in range(1, degree) if degree % d == 0]


class ExtensionDegreeTrait(Trait):
    NAME = "extension_degree"
    DESCRIPTION = (
        "Factorization of the degree $m$ of the field over its prime field. A composite $m$"
        " admits nontrivial intermediate fields; every $m>1$, including prime $m$, also has"
        " the prime field as a proper subfield. These are possible Weil-descent targets, not"
        " by themselves evidence that a GHS transfer is efficient."
    )
    INPUT = {}
    OUTPUT = {
        "degree": (int, "Degree $m$ of the field over its prime field"),
        "factorization": (List[int], "Prime factorization of $m$"),
        "proper_divisor_count": (int, "Number of proper-subfield degrees $d$ with $1 \\le d < m$"),
        "largest_proper_divisor": (
            int,
            "Largest proper-subfield degree, or 0 only when $m=1$",
        ),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Degree of the field extension and its divisor structure.

        A prime field has degree 1 and no proper subfield. A prime-degree binary
        field still contains the prime field, so it reports degree 1 rather than
        incorrectly claiming that no descent target exists. Applicability and
        genus of a concrete descent must be analysed separately.
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
    assert proper_divisors(163) == [1]
    assert proper_divisors(155) == [1, 5, 31]
    assert proper_divisors(176) == [1, 2, 4, 8, 11, 16, 22, 44, 88]


def test_prime_degree_curves_have_only_the_prime_field_target():
    for degree in (113, 163, 233, 283, 409, 571):
        assert proper_divisors(degree) == [1]


def test_composite_degree_curves_do():
    # The X9.62 c2pnb family and the Oakley groups.
    for degree in (155, 176, 185, 208, 272, 304, 368):
        assert proper_divisors(degree) != []
