from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, AnalysisUnavailable, field_value
from dissect.utils.curve_invariants import (
    shioda_invariants, dixmier_ohno_invariants, absolute_ratios,
    higher_binary_invariants, plane_discriminant,
    artin_schreier_orbit,
)


class ModuliInvariantsTrait(DetailedTrait):
    NAME = "moduli_invariants"
    DESCRIPTION = "Igusa-Clebsch, Shioda and Dixmier-Ohno systems with weights and exact absolute ratios."
    INPUT = dict(DetailedTrait.INPUT, max_candidates=(int, "Maximum Artin-Schreier affine orbit candidates"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_candidates=[100000])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "invariant_system": (str, "Named classical invariant system"),
        "normalization": (str, "Construction and scalar normalization"),
        "names": (list, "Invariant coordinate names"),
        "weights": (list, "Degrees in the defining-form coefficients"),
        "values": (list, "Exact field values, in the input field basis"),
        "invariant_count": (int, "Number of weighted coordinates"),
        "vanishing_count": (int, "Number of zero weighted coordinates"),
        "nonzero_chart_index": (int, "Pivot for absolute weight-zero ratios"),
        "absolute_ratios": (list, "Ratios with all numerator and denominator powers specified"),
        "classification_scope": (str, "Geometric moduli; twists need not be distinguished"),
        "complete_classical_system": (int, "Whether this system is complete in the present characteristic"),
        "classification_complete_in_stated_scope": (int, "Completeness in the reported geometric or base-field scope"),
        "canonical_artin_schreier_polynomial": (str, "Canonical reduced polynomial for a single-pole characteristic-two model"),
        "artin_schreier_affine_orbit_size": (int, "Number of reduced affine transforms over the base field"),
        "artin_schreier_affine_stabilizer_order": (int, "Affine stabilizer; each element has two AS lifts"),
    })

    def compute(self, curve, params):
        from sage.all import PolynomialRing
        from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants

        check_budget(curve, params)
        p, g = int(curve.field().characteristic()), curve.genus()
        if curve.model_type() == "CanonicalGenus4":
            from dissect.utils.canonical_genus4 import quadric_matrix

            Q, C = curve.canonical_polynomials()
            names, weights, values = ["quadric_determinant"], [4], [quadric_matrix(Q).det()]
            system, normalization = "Canonical quadric subset", "Determinant of the symmetric quadric matrix; rank and square class are in genus4_canonical_quadric"
            complete = 0
        elif curve.model_type() == "CyclicTrigonal":
            polynomial = curve.cover_polynomial()
            degree = int(polynomial.degree())
            branch_degree = degree+int(degree % 3 != 0)
            names, weights, values = higher_binary_invariants(polynomial, branch_degree)
            system, normalization = "Trigonal branch-divisor subset", "Binary invariants of the unordered branch divisor of the chosen cubic cover"
            complete = 0
        elif curve.is_hyperelliptic():
            if p == 2:
                f, h = curve.hyperelliptic_polynomials()
                orbit = artin_schreier_orbit(f, h, params["max_candidates"])
                polynomial = orbit["polynomial"]
                return {
                    "invariant_system": "Artin-Schreier affine orbit",
                    "normalization": "Remove polynomial AS coboundaries, reduce constants by absolute trace, minimize over all affine x changes",
                    "names": ["a%d" % i for i in range(int(polynomial.degree())+1)],
                    "weights": None, "values": [field_value(c) for c in polynomial],
                    "invariant_count": int(polynomial.degree())+1,
                    "vanishing_count": sum(c == 0 for c in polynomial),
                    "classification_scope": "Base-field isomorphism classes of polynomial single-branch-point Artin-Schreier models, in the fixed input field basis",
                    "complete_classical_system": None, "classification_complete_in_stated_scope": 1,
                    "canonical_artin_schreier_polynomial": str(polynomial),
                    "artin_schreier_affine_orbit_size": orbit["orbit_size"],
                    "artin_schreier_affine_stabilizer_order": orbit["affine_stabilizer_order"],
                }
            f, h = curve.hyperelliptic_polynomials()
            polynomial = f+h**2/4
            if g == 2 and p > 5:
                names, weights = ["I2", "I4", "I6", "I10"], [2, 4, 6, 10]
                values = list(igusa_clebsch_invariants(polynomial))
                system, normalization = "Igusa-Clebsch", "SageMath igusa_clebsch_invariants on f+h^2/4"
                complete = 1
            elif g == 3 and p > 7:
                ring = PolynomialRing(curve.field(), names=("x", "z"))
                x, z = ring.gens()
                binary = sum(coefficient*x**i*z**(8-i) for i, coefficient in enumerate(polynomial))
                values = shioda_invariants(binary)
                names, weights = ["J%d" % i for i in range(2, 11)], list(range(2, 11))
                system, normalization = "Shioda", "Normalized factorial transvectants; J2=(f,f)_8, J3=(f,(f,f)_4)_8; no integral rescaling"
                complete = 1
            else:
                names, weights, values = higher_binary_invariants(polynomial, 2*g+2)
                system, normalization = "Binary-form invariant subset", "Binary discriminant and, when p > 2g+2, normalized transvectants"
                complete = 0
        else:
            degree = int(curve.plane_polynomial().degree())
            if degree == 4 and p > 3:
                names, weights, values = dixmier_ohno_invariants(curve.plane_polynomial())
                system, normalization = "Dixmier-Ohno", "Girard-Kohel section 3 differential normalization; I27=Res(Fx,Fy,Fz)/2^54"
                complete = int(p > 7)
            else:
                names, weights = ["Delta"], [3*(degree-1)**2]
                values = [plane_discriminant(curve.plane_polynomial())]
                system, normalization = "Plane discriminant subset", "Res(Fx,Fy,Fz)/d^(d^2-3d+3)"
                complete = 0
        pivot, ratios = absolute_ratios(values, weights)
        return {
            "invariant_system": system, "normalization": normalization,
            "names": names, "weights": weights, "values": [field_value(v) for v in values],
            "invariant_count": len(values), "vanishing_count": sum(not v for v in values),
            "nonzero_chart_index": pivot, "absolute_ratios": ratios,
            "classification_scope": ("Unordered branch divisor of this trigonal cover, not a complete curve invariant system" if curve.model_type() == "CyclicTrigonal"
                                     else "Canonical quadric only, not a complete curve invariant system" if curve.model_type() == "CanonicalGenus4"
                                     else "Geometric weighted invariants; absolute ratios alone are not a complete isomorphism test and do not classify twists"),
            "complete_classical_system": complete,
            "classification_complete_in_stated_scope": complete,
            "status": "ok" if complete else "partial",
            "reason": None if complete else "The returned coordinates are valid invariants but are not asserted to generate the full invariant ring in this genus and characteristic",
        }
