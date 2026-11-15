"""Validated finite-field models for higher-genus geometric traits.

Subgroup order and cofactor, if supplied, refer to the Jacobian. No elliptic
trace, cardinality, j-invariant or point group law is inferred from them.
"""

from dissect.utils.curve_factory import field_from_dict


def field_element(value, field):
    from sage.all import ZZ
    from dissect.utils.curve_form import dict_to_poly

    if isinstance(value, dict):
        if "poly" in value:
            return dict_to_poly(value["poly"], field)
        value = value["raw"]
        if field.characteristic() == 2 and not field.is_prime_field():
            return field.from_integer(ZZ(value))
    return field(ZZ(value))


class HigherGenusCurve:
    def __init__(self, description):
        from sage.all import Curve, HyperellipticCurve, PolynomialRing, ProjectiveSpace, ZZ

        self._name = description["name"]
        self._category = description["category"]
        self._form = description["form"]
        self._field_type = description["field"]["type"]
        self._field = field_from_dict(description["field"])
        self._standard = description.get("standard", "sim" not in self._name)
        self._example = description.get("example", False)
        self._order = description.get("order")
        self._cofactor = description.get("cofactor")
        if (self._order is None) != (self._cofactor is None):
            raise ValueError("Supply both Jacobian subgroup order and cofactor, or neither")
        if self._order is not None:
            self._order, self._cofactor = ZZ(self._order), ZZ(self._cofactor)
            if self._order <= 0 or self._cofactor <= 0:
                raise ValueError("Jacobian subgroup order and cofactor must be positive")

        params = description["params"]
        if self._form == "Hyperelliptic":
            ring = PolynomialRing(self._field, "x")
            f = ring([field_element(c, self._field) for c in params["f"]])
            h = ring([field_element(c, self._field) for c in params.get("h", [0])])
            self._curve = HyperellipticCurve(f, h, check_squarefree=True)
            self._genus = int(self._curve.genus())
            self._polynomials = f, h
            self._curve_type = "hyperelliptic"
        elif self._form == "SmoothPlane":
            ring = PolynomialRing(self._field, names=("x", "y", "z"))
            polynomial = ring.zero()
            for term in params["polynomial"]:
                powers = tuple(ZZ(p) for p in term["powers"])
                if len(powers) != 3 or any(p < 0 for p in powers):
                    raise ValueError("Plane monomials need three nonnegative powers")
                monomial = field_element(term["coeff"], self._field)
                for variable, power in zip(ring.gens(), powers):
                    monomial *= variable**power
                polynomial += monomial
            if not polynomial or not polynomial.is_homogeneous() or polynomial.degree() < 4:
                raise ValueError("SmoothPlane requires a homogeneous polynomial of degree >= 4")
            ambient = ProjectiveSpace(self._field, 2, names=("x", "y", "z"))
            self._curve = Curve(polynomial, ambient)
            if not self._curve.is_smooth():
                raise ValueError("SmoothPlane requires a geometrically smooth curve")
            degree = int(polynomial.degree())
            self._genus = (degree - 1) * (degree - 2) // 2
            self._curve_type = "non_hyperelliptic"
        else:
            raise ValueError("Unsupported higher-genus form: " + self._form)
        if self._genus < 2:
            raise ValueError("Higher-genus models must have genus >= 2")
        if "genus" in description and ZZ(description["genus"]) != self._genus:
            raise ValueError("Declared genus disagrees with the validated model")

    def name(self):
        return self._name

    def category(self):
        return self._category

    def standard(self):
        return self._standard

    def example(self):
        return self._example

    def field(self):
        return self._field

    def field_type(self):
        return self._field_type

    def q(self):
        return self._field.cardinality()

    def field_bits(self):
        return self.q().nbits()

    def nbits(self):
        return self._order.nbits() if self._order is not None else self.field_bits()

    def order(self):
        return self._order

    def cofactor(self):
        return self._cofactor

    def curve(self):
        return self._curve

    def curve_type(self):
        return self._curve_type

    def genus(self):
        return self._genus

    def is_hyperelliptic(self):
        return self._curve_type == "hyperelliptic"

    def hyperelliptic_polynomials(self):
        if not self.is_hyperelliptic():
            raise ValueError("This model is not hyperelliptic")
        return self._polynomials
