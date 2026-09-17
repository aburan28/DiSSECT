from typing import List

from dissect.traits import Trait

TRAIT_TIMEOUT = 30


def structure_from_invariants(invariants):
    """Summarise the invariant factors of $E(\\mathbb{F}_q)$.

    Sage returns $(n_1, n_2)$ with $n_1 \\mid n_2$ for a bicyclic group and a
    single value for a cyclic one. The first invariant is the interesting part:
    it is 1 exactly when the group is cyclic, and the full 2-torsion (or
    $l$-torsion) being rational is what makes it larger.
    """
    factors = [int(n) for n in invariants if int(n) != 1]
    if not factors:
        return {"invariant_count": 0, "first_invariant": 1, "cyclic": 1}
    if len(factors) == 1:
        return {"invariant_count": 1, "first_invariant": 1, "cyclic": 1}
    return {
        "invariant_count": len(factors),
        "first_invariant": min(factors),
        "cyclic": 0,
    }


class GroupStructureTrait(Trait):
    NAME = "group_structure"
    DESCRIPTION = (
        "Invariant factors of $E(\\mathbb{F}_q)$. A curve whose cofactor is divisible by a square"
        " may be cyclic or $\\mathbb{Z}/n_1 \\times \\mathbb{Z}/n_2$, a structural difference the"
        " cofactor alone cannot express."
    )
    INPUT = {}
    OUTPUT = {
        "invariant_count": (int, "Number of invariant factors greater than 1"),
        "first_invariant": (int, "Smaller invariant factor $n_1$, or 1 when cyclic"),
        "cyclic": (int, "1 when $E(\\mathbb{F}_q)$ is cyclic, 0 otherwise"),
        "invariants": (List[int], "The invariant factors"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Group structure of the curve over its own field.

        Only informative where the cofactor admits it: a prime-order curve is
        cyclic by force, so within a pool of cofactor-1 curves this is constant.
        """
        from dissect.utils.utils import timeout

        invariants = timeout(
            lambda ec: ec.abelian_group().invariants(),
            [curve.ec()],
            timeout_duration=TRAIT_TIMEOUT,
        )
        if isinstance(invariants, str):
            return {
                "invariant_count": None,
                "first_invariant": None,
                "cyclic": None,
                "invariants": None,
            }
        results = structure_from_invariants(invariants)
        results["invariants"] = [int(n) for n in invariants]
        return results


def test_cyclic_group():
    assert structure_from_invariants([7919]) == {
        "invariant_count": 1, "first_invariant": 1, "cyclic": 1
    }


def test_trivial_invariants_are_dropped():
    # Sage may report a trailing 1; it carries no structure.
    assert structure_from_invariants([1, 7919])["cyclic"] == 1
    assert structure_from_invariants([1])["cyclic"] == 1


def test_bicyclic_group():
    result = structure_from_invariants([2, 5000])
    assert result == {"invariant_count": 2, "first_invariant": 2, "cyclic": 0}


def test_first_invariant_is_the_smaller_one():
    assert structure_from_invariants([4, 8])["first_invariant"] == 4


def test_cofactor_four_can_go_either_way():
    # The point of the trait: identical cofactor, different structure.
    cyclic = structure_from_invariants([4 * 1009])
    bicyclic = structure_from_invariants([2, 2 * 1009])
    assert cyclic["cyclic"] == 1 and bicyclic["cyclic"] == 0
