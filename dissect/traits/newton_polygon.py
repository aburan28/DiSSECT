from typing import List
from fractions import Fraction

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, newton_slopes


class NewtonPolygonTrait(DetailedTrait):
    NAME = "newton_polygon"
    DESCRIPTION = "Exact q-normalized Newton slopes, p-rank and isocrystal classifications."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Bound on q^genus"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "vertices": (List[List[int]], "Lower hull vertices with p-adic, unnormalized heights"),
        "slopes": (List[List[int]], "2g exact normalized [numerator, denominator] slopes"),
        "slope_multiplicities": (list, "Distinct slopes with multiplicities"),
        "p_rank": (int, "Multiplicity of slope zero"),
        "ordinary": (int, "1 when p-rank equals genus"),
        "supersingular": (int, "1 when all slopes are 1/2; does not imply superspecial"),
        "slope_symmetry_verified": (int, "Slopes pair to 1 under the principal polarization"),
    })

    def compute(self, curve, params):
        L = l_polynomial(curve, params)
        slopes, vertices = newton_slopes(L, int(curve.field().characteristic()), int(curve.field().degree()))
        return {
            "vertices": [list(p) for p in vertices],
            "slopes": [[s.numerator, s.denominator] for s in slopes],
            "slope_multiplicities": [{"slope": [s.numerator, s.denominator], "multiplicity": slopes.count(s)} for s in sorted(set(slopes))],
            "p_rank": slopes.count(Fraction(0)),
            "ordinary": int(slopes.count(Fraction(0)) == curve.genus()),
            "supersingular": int(all(s == Fraction(1, 2) for s in slopes)),
            "slope_symmetry_verified": int(all(a+b == 1 for a, b in zip(slopes, reversed(slopes)))),
        }
