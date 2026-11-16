from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget
from dissect.utils.cover_analysis import natural_cover, declared_maps, plane_gonality_cover


class CoverRamificationTrait(DetailedTrait):
    NAME = "cover_ramification"
    DESCRIPTION = "Complete local different and branch fibers of natural covers, including wild ramification."
    INPUT = dict(DetailedTrait.INPUT, max_rational_fibers=(int, "Maximum rational target fibers"), max_maps=(int, "Maximum declared maps"), max_map_degree=(int, "Maximum declared-map degree and input exponent"), projection_mode=(int, "1: native; 2: projection through a rational plane point"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_rational_fibers=[64], max_maps=[16], max_map_degree=[256], projection_mode=[1])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "cover_degree": (int, "Degree of this map, not necessarily the gonality"),
        "target_genus": (int, "Genus of the natural-cover target"),
        "base_coordinate": (str, "Separating rational-function coordinate"),
        "extension_coordinate": (str, "Algebraic coordinate in the function-field model"),
        "defining_polynomial": (str, "Equation over the rational function field"),
        "projection_center": (list, "Chosen plane center when projection_mode=2"),
        "projection_coordinate_change": (list, "Columns specifying the new plane basis"),
        "separable": (int, "Separability of the native function-field extension"),
        "different_degree": (int, "Degree of the different, including wild exponents"),
        "riemann_hurwitz_verified": (int, "Checks 2g-2 = -2*degree + different_degree"),
        "branch_closed_points": (int, "Number of closed branch points on the target"),
        "branch_geometric_points": (int, "Geometric branch count, distinct from different degree"),
        "ramified_geometric_points": (int, "Geometric ramified points on the source"),
        "wild_ramified_geometric_points": (int, "Geometric source points with p dividing e"),
        "tame": (int, "1 when every ramification index is prime to p"),
        "branch_fibers": (list, "All branch fibers with e, f, different exponent and local degrees"),
        "rational_fibers": (list, "Base-field rational fibers, including infinity if in the sample"),
        "infinity_fiber": (dict, "Complete fiber above the pole of the separating coordinate"),
        "branch_splitting_extension_degree": (int, "Smallest constant extension making all branch points rational"),
        "ramification_splitting_extension_degree": (int, "Smallest constant extension making all ramified source points rational"),
        "rational_fibers_complete": (int, "Whether all q+1 rational target points were examined"),
        "rational_curve_points_from_fibers": (int, "Independent curve point count when all rational fibers are reported"),
        "discriminant_degree": (int, "Degree of the norm of the different divisor"),
        "geometric_monodromy_group": (str, "C2 for double covers; otherwise not inferred"),
        "galois": (int, "1 for separable double covers; null if not determined"),
        "deck_group_order": (int, "Deck order for double or recognized tame Kummer covers"),
        "galois_constant_extension_degree": (int, "Smallest constant extension making a recognized cover Galois"),
        "declared_maps": (list, "Explicit map equation checks, degrees, target genera and separability"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if params["projection_mode"] not in (1, 2):
            raise ValueError("projection_mode must be 1 or 2")
        projection = None
        if params["projection_mode"] == 2 and curve.model_type() == "SmoothPlane":
            projection = plane_gonality_cover(curve)
        result = natural_cover(curve, params, projection["field"] if projection else None)
        result.update(projection_center=projection["center"] if projection else None,
                      projection_coordinate_change=projection["coordinate_change"] if projection else None,
                      declared_maps=declared_maps(curve, params))
        return result
