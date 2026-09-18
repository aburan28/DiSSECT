from dissect.traits import Trait


def polynomial_profile(exponents):
    """Shape of the polynomial defining an extension field.

    `exponents` are the degrees of its nonzero terms. Three terms is a
    trinomial and five a pentanomial, the two shapes chosen for fast reduction;
    anything longer is a different kind of object, in practice the polynomial
    attached to a normal basis rather than a polynomial basis. The second
    highest degree is what governs reduction cost, since it sets how far a
    carry out of the top propagates.
    """
    degrees = sorted((int(e) for e in exponents), reverse=True)
    return {
        "terms": len(degrees),
        "reduction_degree": degrees[1] if len(degrees) > 1 else 0,
    }


class FieldRepresentationTrait(Trait):
    NAME = "field_representation"
    DESCRIPTION = (
        "Shape of the polynomial defining an extension field: how many terms it has, and the"
        " second highest degree, which governs reduction cost. Trinomials and pentanomials are"
        " the fast polynomial-basis choices; a longer polynomial marks a normal basis instead."
    )
    INPUT = {}
    OUTPUT = {
        "terms": (int, "Number of nonzero terms in the defining polynomial"),
        "reduction_degree": (int, "Second highest degree present"),
    }
    DEFAULT_PARAMS = {}

    def compute(self, curve, params):
        """Defining-polynomial shape, for extension fields only.

        A prime field has no such polynomial to describe. `prime_shape` is the
        counterpart there; it is silent on binary fields for the mirror-image
        reason, since $q = 2^m$ is a single bit however the field is built.
        """
        if curve.is_over_prime():
            return {"terms": None, "reduction_degree": None}
        return polynomial_profile(curve.field().modulus().exponents())


def test_trinomial():
    # sect233k1: x^233 + x^74 + 1
    assert polynomial_profile([233, 74, 0]) == {"terms": 3, "reduction_degree": 74}


def test_pentanomial():
    # B-163 and K-163: x^163 + x^7 + x^6 + x^3 + 1
    assert polynomial_profile([163, 7, 6, 3, 0]) == {"terms": 5, "reduction_degree": 7}


def test_order_of_exponents_does_not_matter():
    assert polynomial_profile([0, 74, 233]) == polynomial_profile([233, 74, 0])


def test_normal_basis_polynomials_are_long():
    # c2onb191v4 is given by a 15-term polynomial, the signature of a normal
    # basis rather than the trinomial or pentanomial of a polynomial basis.
    onb = [191, 190, 188, 184, 176, 160, 128, 64, 32, 16, 8, 4, 2, 1, 0]
    profile = polynomial_profile(onb)
    assert profile["terms"] == 15
    assert profile["terms"] not in (3, 5)


def test_degenerate_input():
    assert polynomial_profile([7]) == {"terms": 1, "reduction_degree": 0}
