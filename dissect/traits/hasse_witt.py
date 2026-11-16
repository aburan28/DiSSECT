from typing import List

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, function_field, field_matrix


class HasseWittTrait(DetailedTrait):
    NAME = "hasse_witt"
    DESCRIPTION = "Cartier and dual Hasse-Witt matrices, stable semilinear p-rank and a-number."
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "cartier_matrix": (list, "Column images in the reported regular differential basis"),
        "hasse_witt_matrix": (list, "Frobenius-semilinear Serre-dual matrix"),
        "differential_basis": (List[str], "Basis used for the matrices"),
        "cartier_rank": (int, "Rank of one Cartier application"),
        "stable_cartier_ranks": (List[int], "Ranks of semilinear iterates 1 through genus"),
        "p_rank": (int, "Stable rank, not necessarily the first matrix rank"),
        "a_number": (int, "genus minus the Cartier rank"),
        "ordinary": (int, "p-rank equals genus"),
        "superspecial": (int, "Cartier zero: Jacobian geometrically superspecial"),
        "nilpotence_index": (int, "First zero Cartier iterate; null when not nilpotent"),
    })

    def compute(self, curve, params):
        from sage.all import matrix, identity_matrix

        check_budget(curve, params)
        field = curve.field()
        p, e, g = int(field.characteristic()), int(field.degree()), curve.genus()
        inverse_frobenius = p**(e-1)
        if curve.is_hyperelliptic() and p != 2:
            from dissect.utils.curve_analysis import AnalysisUnavailable

            f, h = curve.hyperelliptic_polynomials()
            if int((f+h**2/4).degree())*((p-1)//2) > 2000000:
                raise AnalysisUnavailable("Cartier coefficient polynomial exceeds the term budget", "resource_limit")
            coefficients = (f+h**2/4)**((p-1)//2)
            C = matrix(field, g, lambda i, j: coefficients[p*(i+1)-(j+1)]**inverse_frobenius if p*(i+1) >= j+1 else 0)
            basis = ["x^%d dx/(y+h(x)/2)" % i for i in range(g)]
        else:
            L = function_field(curve)
            V, to_differential, to_vector = L.space_of_holomorphic_differentials()
            basis = [to_differential(v) for v in V.basis()]
            C = matrix(field, [to_vector(w.cartier()) for w in basis]).transpose()
        iterate = identity_matrix(field, g)
        ranks = []
        for _ in range(g):
            iterate = C * iterate.apply_map(lambda a: a**inverse_frobenius)
            ranks.append(int(iterate.rank()))
        H = C.apply_map(lambda a: a**p).transpose()
        return {
            "cartier_matrix": field_matrix(C), "hasse_witt_matrix": field_matrix(H),
            "differential_basis": [str(w) for w in basis],
            "cartier_rank": int(C.rank()), "stable_cartier_ranks": ranks,
            "p_rank": ranks[-1], "a_number": g-int(C.rank()),
            "ordinary": int(ranks[-1] == g), "superspecial": int(C.is_zero()),
            "nilpotence_index": ranks.index(0)+1 if 0 in ranks else None,
        }
