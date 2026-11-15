"""Dispatch database models without treating a Jacobian as an elliptic curve."""


HIGHER_GENUS_FORMS = ("Hyperelliptic", "SmoothPlane")


def field_from_dict(description):
    from sage.all import GF, ZZ
    from dissect.utils.curve_form import dict_to_poly

    if description["type"] == "Prime":
        return GF(ZZ(description["p"]), proof=False)
    base = ZZ(2) if description["type"] == "Binary" else ZZ(description["base"])
    degree = ZZ(description["degree"])
    modulus = dict_to_poly(description["poly"], GF(base)["w"])
    return GF(base**degree, "w", modulus, proof=False)


def curve_from_dict(description):
    if description["form"] in HIGHER_GENUS_FORMS:
        from dissect.utils.higher_genus_curve import HigherGenusCurve

        return HigherGenusCurve(description)
    from dissect.utils.custom_curve import CustomCurve

    return CustomCurve(description)
