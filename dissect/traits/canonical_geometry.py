from dissect.traits import Trait


def canonical_profile(genus, hyperelliptic):
    """Geometry over the algebraic closure; hyperelliptic means genus >= 2."""
    if genus < 1 or (hyperelliptic and genus < 2) or (genus == 2 and not hyperelliptic):
        raise ValueError("Inconsistent genus and curve family")
    elliptic = genus == 1
    return {
        "genus": genus,
        "holomorphic_differentials": genus,
        "jacobian_dimension": genus,
        "jacobian_is_curve": int(elliptic),
        "canonical_projective_dimension": genus - 1,
        "canonical_map_degree": None if elliptic else (2 if hyperelliptic else 1),
        "canonical_image_degree": None if elliptic else (genus - 1 if hyperelliptic else 2 * genus - 2),
        "canonical_embedding": int(not elliptic and not hyperelliptic),
        "degree_two_map": int(elliptic or hyperelliptic),
        "unique_degree_two_map": int(hyperelliptic),
        "plane_quartic": int(genus == 3 and not hyperelliptic),
    }


class CanonicalGeometryTrait(Trait):
    NAME = "canonical_geometry"
    DESCRIPTION = "Genus, differentials, Jacobian dimension and canonical-map geometry."
    SUPPORTED_CURVE_TYPES = ("elliptic", "hyperelliptic", "non_hyperelliptic")
    OUTPUT = {
        "genus": (int, "Geometric genus"),
        "holomorphic_differentials": (int, "Dimension of the regular differential space"),
        "jacobian_dimension": (int, "Dimension of the Jacobian"),
        "jacobian_is_curve": (int, "1 when the curve is its Jacobian after choosing an origin"),
        "canonical_projective_dimension": (int, "Target projective dimension of the canonical map"),
        "canonical_map_degree": (int, "Generic degree; null for the constant genus-one map"),
        "canonical_image_degree": (int, "Degree of the canonical image curve; null in genus one"),
        "canonical_embedding": (int, "1 when the canonical map is an embedding"),
        "degree_two_map": (int, "1 when a degree-two map to the projective line exists geometrically"),
        "unique_degree_two_map": (int, "1 when that map is unique up to target automorphisms"),
        "plane_quartic": (int, "1 when the canonical model is a smooth plane quartic"),
    }

    def compute(self, curve, params):
        return canonical_profile(int(curve.genus()), curve.is_hyperelliptic())
