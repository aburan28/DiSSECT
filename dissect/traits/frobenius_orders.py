from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, AnalysisUnavailable


class FrobeniusOrdersTrait(DetailedTrait):
    NAME = "frobenius_orders"
    DESCRIPTION = "CM fields and Frobenius/Verschiebung orders for simple ordinary Jacobians."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Point-count budget"), max_order_degree=(int, "Number-field degree bound"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], max_order_degree=[12])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "frobenius_polynomial": (list, "Monic irreducible polynomial, coefficients ascending"),
        "cm_field_degree": (int, "Degree of the commutative base-field endomorphism algebra"),
        "cm_field_discriminant": (int, "Discriminant of the maximal order of Q(pi)"),
        "cm_field_signature": (list, "Real embeddings and complex pairs"),
        "frobenius_order_discriminant": (int, "Discriminant of Z[pi]"),
        "frobenius_order_index": (int, "Index of Z[pi] in the maximal order"),
        "weil_order_discriminant": (int, "Discriminant of Z[pi,q/pi]"),
        "weil_order_index": (int, "Index of Z[pi,q/pi] in the maximal order"),
        "verschiebung": (str, "q/pi in the reported number-field basis"),
        "real_subfield_polynomial": (list, "Minimal polynomial of pi+q/pi"),
        "real_subfield_discriminant": (int, "Maximal real-subfield discriminant"),
        "complex_conjugation_verified": (int, "pi -> q/pi is an involution fixing pi+q/pi"),
        "discriminant_index_identities_verified": (int, "Both order discriminants equal field discriminant times index squared"),
        "actual_endomorphism_order_determined": (int, "0 unless the containing order bounds coincide"),
        "actual_endomorphism_order_maximal": (int, "1 when the Weil order is already maximal; otherwise null"),
        "endomorphism_order_scope": (str, "Base-field endomorphisms; no unsupported geometric-ring claim"),
    })

    def compute(self, curve, params):
        from sage.all import PolynomialRing, QQ, NumberField, ZZ

        coefficients = [int(c) for c in l_polynomial(curve, params)]
        if coefficients[curve.genus()] % int(curve.field().characteristic()) == 0:
            raise AnalysisUnavailable("Requires an ordinary Jacobian", "not_applicable")
        ring = PolynomialRing(QQ, "T")
        polynomial = ring(list(reversed(coefficients)))
        if not polynomial.is_irreducible():
            raise AnalysisUnavailable("Requires a simple ordinary Jacobian with irreducible Frobenius polynomial", "not_applicable")
        if polynomial.degree() > params["max_order_degree"]:
            raise AnalysisUnavailable("Number-field degree exceeds max_order_degree", "resource_limit")
        field = NumberField(polynomial, "pi")
        pi, q = field.gen(), ZZ(curve.q())
        verschiebung = q/pi
        conjugation = field.hom([verschiebung], field)
        if conjugation(verschiebung) != pi or conjugation(pi+verschiebung) != pi+verschiebung:
            raise ValueError("The Weil-field conjugation failed its involution identities")
        if field.signature() != (0, curve.genus()):
            raise ValueError("A simple ordinary Weil field must be a CM field")
        frobenius_order = field.order([pi])
        weil_order = field.order([pi, verschiebung])
        maximal = field.maximal_order()
        field_discriminant = ZZ(field.discriminant())
        indices = [ZZ(order.index_in(maximal)) for order in (frobenius_order, weil_order)]
        for order, index in zip((frobenius_order, weil_order), indices):
            if order.discriminant() != field_discriminant*index*index:
                raise ValueError("A CM-order discriminant/index identity failed")
        real_polynomial = (pi+verschiebung).minpoly()
        if real_polynomial.degree() != curve.genus():
            raise ValueError("The real Weil subfield has unexpected degree")
        real_field = NumberField(real_polynomial, "u")
        return {"frobenius_polynomial": [int(c) for c in polynomial], "cm_field_degree": int(field.degree()),
                "cm_field_discriminant": int(field_discriminant), "cm_field_signature": list(field.signature()),
                "frobenius_order_discriminant": int(frobenius_order.discriminant()), "frobenius_order_index": int(indices[0]),
                "weil_order_discriminant": int(weil_order.discriminant()), "weil_order_index": int(indices[1]),
                "verschiebung": str(verschiebung), "real_subfield_polynomial": [int(c) for c in real_polynomial],
                "real_subfield_discriminant": int(real_field.discriminant()), "complex_conjugation_verified": 1,
                "discriminant_index_identities_verified": 1, "actual_endomorphism_order_determined": int(indices[1] == 1),
                "actual_endomorphism_order_maximal": 1 if indices[1] == 1 else None,
                "endomorphism_order_scope": "Z[pi,q/pi] is contained in End_Fq(J), which is contained in the maximal CM order; intermediate order is not inferred"}
