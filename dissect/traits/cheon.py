"""Conditional Brown--Gallant--Cheon attack screening.

This is deliberately a protocol-level trait. A divisor of the prime subgroup
order ``r - 1`` or ``r + 1`` does not speed up an ordinary, single-instance
ECDLP. It matters when an attacker receives suitable auxiliary powers or can
query a service which repeatedly applies one static secret scalar.
"""

import math
from typing import List

from dissect.traits import Trait


TRAIT_TIMEOUT = 30
MAX_DIVISORS = 1_000_000


def _log2_integer(value):
    """Accurate log2 for positive integers without converting them to float."""
    if value <= 0:
        raise ValueError("log2 is defined here only for positive integers")
    exponent = value.bit_length() - 1
    return exponent + math.log2(value / (1 << exponent))


def _log2_sum(left, right):
    """Return log2(2**left + 2**right) without overflow."""
    high, low = max(left, right), min(left, right)
    return high + math.log2(1.0 + 2.0 ** (low - high))


def divisor_count(factorization):
    """Number of positive divisors represented by ``[(prime, exponent), ...]``."""
    count = 1
    for _, exponent in factorization:
        count *= int(exponent) + 1
    return count


def divisors_from_factorization(factorization, limit=MAX_DIVISORS):
    """Enumerate divisors exactly, refusing unexpectedly large searches."""
    count = divisor_count(factorization)
    if count > limit:
        raise ValueError(f"factorization has {count} divisors (limit {limit})")

    divisors = [1]
    for prime, exponent in factorization:
        prime, exponent = int(prime), int(exponent)
        previous = tuple(divisors)
        power = 1
        for _ in range(exponent):
            power *= prime
            divisors.extend(divisor * power for divisor in previous)
    return divisors


def oracle_attack_profile(order, divisor):
    """Cost model for the conditional static-scalar oracle attack.

    The simplified Brown--Gallant--Cheon tradeoff uses ``u`` oracle queries and
    about ``sqrt(r/u)`` offline group operations for a divisor ``u`` of
    ``r - 1`` or ``r + 1``. The returned total is their sum, not merely the
    larger phase.
    """
    order, divisor = int(order), int(divisor)
    if order <= 3 or not 1 < divisor < order:
        raise ValueError("require 1 < divisor < subgroup order")

    query_bits = _log2_integer(divisor)
    offline_bits = 0.5 * (_log2_integer(order) - query_bits)
    total_bits = _log2_sum(query_bits, offline_bits)
    return {
        "divisor": divisor,
        "query_bits": query_bits,
        "offline_bits": offline_bits,
        "total_bits": total_bits,
    }


def best_oracle_attack(
    order, factorization, limit=MAX_DIVISORS, max_query_bits=None
):
    """Return the lowest-cost profile within an optional oracle-query budget."""
    if max_query_bits is not None and max_query_bits < 1:
        raise ValueError("max_query_bits must be at least 1 or None")
    if max_query_bits is None or max_query_bits >= int(order).bit_length():
        max_divisor = None
    else:
        max_divisor = 1 << max_query_bits
    profiles = (
        oracle_attack_profile(order, divisor)
        for divisor in divisors_from_factorization(factorization, limit)
        if 1 < divisor < order
        and (max_divisor is None or divisor <= max_divisor)
    )
    return min(profiles, key=lambda profile: profile["total_bits"])


class CheonTrait(Trait):
    NAME = "cheon"
    DESCRIPTION = (
        "Conditional Brown--Gallant--Cheon screen for factors of the prime-subgroup"
        " order $r\\pm1$. It reports the best tradeoff between $u$ static-scalar"
        " oracle queries and about $\\sqrt{r/u}$ offline operations. This is not a"
        " generic ECDLP security estimate: without the auxiliary-input/oracle"
        " precondition, Pollard rho remains the applicable comparison."
    )
    INPUT = {
        "max_query_bits": (
            int,
            "Maximum $\\log_2$ oracle queries; -1 searches the unconstrained optimum",
        )
    }
    OUTPUT = {
        "subgroup_order_prime": (int, "1 when the supplied subgroup order is prime"),
        "r_minus_1_factorization": (List[int], "Factorization of $r-1$"),
        "r_plus_1_factorization": (List[int], "Factorization of $r+1$"),
        "factorization_complete": (int, "1 when both neighboring orders were factored"),
        "search_complete": (int, "1 when every divisor of each completed factorization was searched"),
        "best_sign": (str, "Neighbor supplying the best known divisor: '-' or '+'"),
        "best_divisor": (int, "Best known oracle-query divisor $u$"),
        "query_bits": (float, "$\\log_2 u$ oracle queries"),
        "offline_bits": (float, "$\\log_2\\sqrt{r/u}$ offline group operations"),
        "total_bits": (float, "$\\log_2(u+\\sqrt{r/u})$ conditional attack cost"),
        "rho_bits": (float, "$\\log_2\\sqrt r$ generic baseline, ignoring constants"),
        "advantage_bits": (float, "Generic baseline minus conditional total cost"),
    }
    DEFAULT_PARAMS = {"max_query_bits": [-1]}

    def compute(self, curve, params):
        from sage.all import ZZ
        from dissect.utils.utils import Factorization

        order = int(curve.order())
        max_query_bits = int(params["max_query_bits"])
        if max_query_bits == 0 or max_query_bits < -1:
            raise ValueError("max_query_bits must be -1 (unconstrained) or at least 1")
        query_budget = None if max_query_bits == -1 else max_query_bits
        prime = bool(ZZ(order).is_prime(proof=False))
        result = {
            "subgroup_order_prime": int(prime),
            "r_minus_1_factorization": None,
            "r_plus_1_factorization": None,
            "factorization_complete": 0,
            "search_complete": 0,
            "best_sign": None,
            "best_divisor": None,
            "query_bits": None,
            "offline_bits": None,
            "total_bits": None,
            "rho_bits": 0.5 * _log2_integer(order) if order > 0 else None,
            "advantage_bits": None,
        }
        if not prime:
            return result

        candidates = []
        complete = True
        search_complete = True
        for sign, neighbor, output in (
            ("-", order - 1, "r_minus_1_factorization"),
            ("+", order + 1, "r_plus_1_factorization"),
        ):
            factors = Factorization(neighbor, timeout_duration=TRAIT_TIMEOUT)
            result[output] = factors.factorization()
            if factors.timeout():
                complete = False
                continue
            try:
                profile = best_oracle_attack(
                    order,
                    factors.factorization(unpack=False),
                    max_query_bits=query_budget,
                )
            except ValueError:
                search_complete = False
                continue
            profile["sign"] = sign
            candidates.append(profile)

        result["factorization_complete"] = int(complete)
        result["search_complete"] = int(search_complete and complete)
        if not candidates:
            return result

        best = min(candidates, key=lambda profile: profile["total_bits"])
        result.update(
            {
                "best_sign": best["sign"],
                "best_divisor": best["divisor"],
                "query_bits": best["query_bits"],
                "offline_bits": best["offline_bits"],
                "total_bits": best["total_bits"],
                "advantage_bits": result["rho_bits"] - best["total_bits"],
            }
        )
        return result


def test_oracle_attack_profile_balances_the_two_phases():
    profile = oracle_attack_profile(2**90 + 1, 2**30)
    assert math.isclose(profile["query_bits"], 30.0)
    assert math.isclose(profile["offline_bits"], 30.0, abs_tol=1e-12)
    assert math.isclose(profile["total_bits"], 31.0, abs_tol=1e-12)


def test_best_oracle_attack_checks_composite_divisors():
    # The best divisor is composite, not an individual prime factor.
    result = best_oracle_attack(1009, [(2, 4), (3, 2), (7, 1)])
    assert result["divisor"] == 6


def test_best_oracle_attack_ignores_r_plus_1_itself():
    result = best_oracle_attack(101, [(2, 1), (3, 1), (17, 1)])
    assert result["divisor"] == 3


def test_best_oracle_attack_honors_query_budget():
    result = best_oracle_attack(
        1009, [(2, 4), (3, 2), (7, 1)], max_query_bits=2
    )
    assert result["divisor"] == 4


def test_divisor_search_is_bounded():
    try:
        divisors_from_factorization([(2, 20)], limit=20)
    except ValueError as error:
        assert "21 divisors" in str(error)
    else:
        raise AssertionError("oversized divisor search was not rejected")
