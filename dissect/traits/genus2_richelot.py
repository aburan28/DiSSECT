from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable
from dissect.utils.richelot import richelot_targets


class Genus2RichelotTrait(DetailedTrait):
    NAME = "genus2_richelot"
    DESCRIPTION = "All fifteen maximal isotropic two-torsion kernels, Richelot targets and split degenerations."
    GENUS = 2
    SUPPORTED_CURVE_TYPES = ("hyperelliptic",)
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Source point-count budget"), max_target_point_work=(int, "Per-target q_search squared verification bound"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], max_target_point_work=[128])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "geometric_kernel_count": (int, "Fifteen maximal isotropic subgroups in geometric J[2]"),
        "base_field_kernel_count": (int, "Frobenius-stable kernels, even when individual pairs are permuted"),
        "kernel_order": (int, "Four elements in a (2,2)-kernel"), "polarized_isogeny_degree": (int, "Degree four as an abelian-surface isogeny"),
        "kernel_maximally_isotropic": (int, "Pair partitions give Lagrangians for the Weil pairing"),
        "branch_splitting_extension_degree": (int, "Constant extension used for explicit equations"),
        "branch_points": (list, "Homogeneous branch coordinates, including infinity"),
        "branch_frobenius": (dict, "Frobenius on branch points"), "kernel_frobenius": (dict, "Frobenius on the fifteen pair partitions"),
        "target_field": (dict, "Field and embedding for explicit targets"),
        "kernel_partitions": (list, "Quadratic factors, determinants, bracket certificates, and targets"),
        "genus_two_target_count": (int, "Nondegenerate Richelot quotient Jacobians"),
        "split_product_kernel_count": (int, "Determinant-zero quotients with product principal polarization"),
        "distinct_target_equation_count": (int, "Literal equations, not isomorphism classes"),
        "targets_complete": (int, "All kernels classified and every nondegenerate target equation constructed"),
        "point_count_verified_targets": (int, "Targets additionally checked by independent point counts"),
        "target_descent_scope": (str, "Explicit field versus minimal kernel field; split components not constructed"),
    })

    def compute(self, curve, params):
        check_budget(curve, params)
        if curve.field().characteristic() == 2:
            raise AnalysisUnavailable("Classical separable Richelot two-isogenies require odd characteristic")
        return richelot_targets(curve, params)
