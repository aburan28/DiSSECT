from dissect.traits import Trait


def relative_size(discriminant, q):
    """Bit length of $|D|$ against the bit length of the field.

    A curve built by complex multiplication starts from a chosen, tiny
    discriminant and solves for a field, so the ratio is near zero. A curve
    found by hashing a seed takes whatever discriminant it gets, and that is
    the size of the field.
    """
    return abs(discriminant).bit_length() / q.bit_length()


class CMDiscriminantSizeTrait(Trait):
    NAME = "cm_discriminant_size"
    DESCRIPTION = (
        "Size of the CM discriminant relative to the field. `discriminant` already reports $D$"
        " itself, but as a several-hundred-bit signed integer it is hard to compare; the ratio"
        " of bit lengths separates CM-constructed curves, which is near 0, from seed-generated"
        " ones, which is near 1."
    )
    INPUT = {}
    OUTPUT = {
        "discriminant_bitlen": (int, "Bit length of $|D|$"),
        "relative_size": (float, "Bit length of $|D|$ over bit length of $q$"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Relative size of the CM discriminant.

        `cm_discriminant` factors the Frobenius discriminant and can time out,
        in which case it hands back a message rather than a number.
        """
        discriminant = curve.cm_discriminant()
        if not isinstance(discriminant, (int, type(curve.q()))):
            return {"discriminant_bitlen": None, "relative_size": None}

        discriminant, q = int(discriminant), int(curve.q())
        return {
            "discriminant_bitlen": abs(discriminant).bit_length(),
            "relative_size": relative_size(discriminant, q),
        }


def test_cm_constructed_curves_are_near_zero():
    # BN and BLS curves are built at j = 0, so D = -3 whatever the field size.
    assert relative_size(-3, 2**158) < 0.02
    assert relative_size(-3, 2**381) < 0.01
    # MNT curves start from a small chosen discriminant.
    assert relative_size(-19, 2**170) < 0.04
    assert relative_size(-91, 2**159) < 0.05


def test_seed_generated_curves_are_near_one():
    # A discriminant of the same size as the field, which is what hashing gives.
    d = -(2**255 + 12345)
    assert 0.95 < relative_size(d, 2**256) < 1.05


def test_sign_is_ignored():
    # The CM discriminant is negative for an ordinary curve; only size matters.
    assert relative_size(-3, 2**158) == relative_size(3, 2**158)


def test_separation_is_wide():
    # The point of the trait: two orders of magnitude between the regimes.
    cm = relative_size(-3, 2**381)
    seeded = relative_size(-(2**380), 2**381)
    assert seeded / cm > 100
