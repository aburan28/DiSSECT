import math

from dissect.traits import Trait


def automorphism_group_order(characteristic, field_size, j_invariant):
    """Order of the automorphism group of a curve with the given $j$-invariant.

    Away from characteristic 2 and 3 the extra automorphisms at $j = 0$ and
    $j = 1728$ are only defined over the base field when it contains the
    corresponding roots of unity, so the field size decides. In characteristic
    2 and 3 the $j = 0$ curves are supersingular and the counts below are those
    over the algebraic closure.
    """
    if characteristic == 2:
        return 24 if j_invariant == 0 else 2
    if characteristic == 3:
        return 12 if j_invariant == 0 else 2
    if j_invariant == 0:
        return 6 if field_size % 3 == 1 else 2
    if j_invariant == 1728 % characteristic:
        return 4 if field_size % 4 == 1 else 2
    return 2


def rho_speedup_bits(group_order):
    """Margin lost relative to a generic curve, whose only automorphisms are $\\pm 1$."""
    return 0.5 * math.log2(group_order / 2)


class AutomorphismsTrait(Trait):
    NAME = "automorphisms"
    DESCRIPTION = (
        "Order of the automorphism group. Beyond the generic $\\pm 1$, the extra automorphisms"
        " at $j = 0$ and $j = 1728$ enlarge the Pollard rho equivalence classes and so shave"
        " $\\log_2\\sqrt{|\\mathrm{Aut}|/2}$ bits off the search."
    )
    INPUT = {}
    OUTPUT = {
        "automorphism_group_order": (int, "Order of $\\mathrm{Aut}(E)$"),
        "rho_speedup_bits": (float, "$\\log_2\\sqrt{|\\mathrm{Aut}|/2}$"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Automorphism count derived from the $j$-invariant and the field."""
        from sage.all import ZZ

        field = curve.field()
        order = automorphism_group_order(
            int(field.characteristic()),
            int(curve.q()),
            int(ZZ(curve.j_invariant())),
        )
        return {
            "automorphism_group_order": order,
            "rho_speedup_bits": rho_speedup_bits(order),
        }


def test_generic_curve():
    # A 256-bit prime field, j neither 0 nor 1728.
    assert automorphism_group_order(23, 23, 5) == 2


def test_j_zero_needs_cube_roots_of_unity():
    # secp256k1 has j = 0 over a field with q = 1 mod 3, giving the GLV endomorphism.
    secp256k1_p = 2**256 - 2**32 - 977
    assert secp256k1_p % 3 == 1
    assert automorphism_group_order(secp256k1_p, secp256k1_p, 0) == 6
    # Without cube roots of unity the extra automorphisms are not defined over the field.
    assert automorphism_group_order(5, 5, 0) == 2


def test_j_1728_needs_fourth_roots_of_unity():
    assert automorphism_group_order(13, 13, 1728 % 13) == 4
    assert automorphism_group_order(7, 7, 1728 % 7) == 2


def test_small_characteristic():
    assert automorphism_group_order(2, 2**163, 0) == 24
    assert automorphism_group_order(3, 3**97, 0) == 12
    # The Koblitz curves are ordinary with j != 0, so their speedup is Frobenius, not Aut.
    assert automorphism_group_order(2, 2**163, 1) == 2


def test_rho_speedup_bits():
    assert rho_speedup_bits(2) == 0.0
    assert math.isclose(rho_speedup_bits(6), 0.7925, abs_tol=1e-4)
    assert math.isclose(rho_speedup_bits(4), 0.5, abs_tol=1e-9)
