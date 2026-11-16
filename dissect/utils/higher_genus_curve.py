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
        from sage.all import Curve, HyperellipticCurve, PolynomialRing, ProjectiveSpace, FunctionField, ZZ

        self._name = description["name"]
        self._category = description["category"]
        self._form = description["form"]
        self._field_type = description["field"]["type"]
        self._field_description = description["field"]
        self._field = field_from_dict(description["field"])
        self._standard = description.get("standard", "sim" not in self._name)
        self._example = description.get("example", False)
        self._analysis = description.get("analysis", {})
        if not isinstance(self._analysis, dict) or not isinstance(self._analysis.get("maps", []), list):
            raise ValueError("analysis must be an object containing a maps list")
        if any(not isinstance(mapping, dict) for mapping in self._analysis.get("maps", [])):
            raise ValueError("Each declared map must be an object")
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
        elif self._form == "CyclicTrigonal":
            if self._field.characteristic() == 3:
                raise ValueError("CyclicTrigonal requires characteristic different from three")
            ring = PolynomialRing(self._field, "x")
            f = ring([field_element(c, self._field) for c in params["f"]])
            if f.degree() < 4 or not f.is_squarefree():
                raise ValueError("CyclicTrigonal requires a squarefree polynomial of degree >= 4")
            K = FunctionField(self._field, "x")
            equation_ring = PolynomialRing(K, "y")
            L = K.extension(equation_ring.gen()**3-f(K.gen()), "y")
            self._curve = L
            self._cover_function_field = L
            self._cover_polynomial = f
            self._genus = int(L.genus())
            self._curve_type = "non_hyperelliptic"
        elif self._form == "CanonicalGenus4":
            from dissect.utils.canonical_genus4 import homogeneous_form, quadric_matrix

            if self._field.characteristic() == 2:
                raise ValueError("CanonicalGenus4 currently requires odd characteristic")
            ring = PolynomialRing(self._field, names=("X0", "X1", "X2", "X3"))
            Q = homogeneous_form(ring, params["quadric"], 2)
            C = homogeneous_form(ring, params["cubic"], 3)
            if quadric_matrix(Q).rank() not in (3, 4):
                raise ValueError("CanonicalGenus4 requires a quadric of rank three or four")
            self._curve = Curve([Q, C], ProjectiveSpace(self._field, 3, names=ring.variable_names()))
            if self._curve.dimension() != 1 or not self._curve.is_smooth():
                raise ValueError("CanonicalGenus4 requires a geometrically smooth quadric-cubic intersection")
            self._canonical_polynomials = Q, C
            self._genus = 4
            self._curve_type = "non_hyperelliptic"
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
            self._plane_polynomial = polynomial
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

    def field_description(self):
        return self._field_description

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

    def model_type(self):
        return self._form

    def genus(self):
        return self._genus

    def is_hyperelliptic(self):
        return self._curve_type == "hyperelliptic"

    def hyperelliptic_polynomials(self):
        if not self.is_hyperelliptic():
            raise ValueError("This model is not hyperelliptic")
        return self._polynomials

    def plane_polynomial(self):
        if self._form != "SmoothPlane":
            raise ValueError("This model is not a smooth plane curve")
        return self._plane_polynomial

    def analysis(self):
        return self._analysis

    def cover_function_field(self):
        return self._cover_function_field

    def cover_polynomial(self):
        return self._cover_polynomial

    def canonical_polynomials(self):
        return self._canonical_polynomials
