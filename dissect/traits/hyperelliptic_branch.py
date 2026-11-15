from typing import List

from dissect.traits import Trait


class HyperellipticBranchTrait(Trait):
    NAME = "hyperelliptic_branch"
    DESCRIPTION = "Tame branch divisor and Frobenius orbits of a hyperelliptic double cover."
    SUPPORTED_CURVE_TYPES = ("hyperelliptic",)
    OUTPUT = {
        "geometric_branch_points": (int, "Number of geometric branch points"),
        "rational_branch_points": (int, "Number of base-field rational branch points"),
        "branch_at_infinity": (int, "1 when infinity is a branch point"),
        "branch_closed_points": (int, "Number of Frobenius orbits in the branch divisor"),
        "branch_orbit_degrees": (List[int], "Sorted Frobenius orbit lengths, including infinity"),
    }

    def compute(self, curve, params):
        if curve.field().characteristic() == 2:
            return dict.fromkeys(self.OUTPUT)
        f, h = curve.hyperelliptic_polynomials()
        branch_polynomial = h**2 + 4 * f
        degrees = [int(factor.degree()) for factor, _ in branch_polynomial.factor()]
        infinity = int(branch_polynomial.degree() % 2)
        if infinity:
            degrees.append(1)
        return {
            "geometric_branch_points": sum(degrees),
            "rational_branch_points": degrees.count(1),
            "branch_at_infinity": infinity,
            "branch_closed_points": len(degrees),
            "branch_orbit_degrees": sorted(degrees),
        }
