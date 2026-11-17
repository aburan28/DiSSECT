from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, field_matrix, AnalysisUnavailable
from dissect.utils.jacobian_profiles import torsion_constraints, hyperelliptic_two_torsion


class JacobianTorsionTrait(DetailedTrait):
    NAME = "jacobian_torsion"
    DESCRIPTION = "Prime-to-p Jacobian torsion ranks and extension bounds, with exact hyperelliptic two-torsion."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Point-count budget"), ell=(int, "Prime distinct from the characteristic"),
                 extension_bound=(int, "Largest tested constant extension"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], ell=[2], extension_bound=[12])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "ell": (int, "Torsion prime"), "module_dimension": (int, "Geometric prime-to-p torsion dimension 2g"),
        "frobenius_factors_mod_ell": (list, "Irreducible primary factors with multiplicities"),
        "rank_information_source": (str, "Exact even-branch-subset module or characteristic-polynomial bounds"),
        "frobenius_matrix": (list, "Exact branch-subset quotient matrix when available"),
        "extension_ranks": (list, "Rational ranks or rigorous Jordan-form bounds in each extension"),
        "base_rational_rank": (int, "Exact rational rank if determined"),
        "base_rank_lower_bound": (int, "Guaranteed rational torsion rank"),
        "base_rank_upper_bound": (int, "Largest rank allowed by the known module data"),
        "least_nontrivial_extension_degree": (int, "Exact first nonzero extension if found or known from the module"),
        "semisimple_frobenius_order": (int, "Order after removing the unipotent part, when determined"),
        "full_torsion_extension_degree": (int, "Exact full-torsion extension only when certified"),
        "full_torsion_extension_upper_bound": (int, "A certified degree containing the full torsion"),
        "unipotent_order_upper_bound": (int, "Largest possible ell-power part from primary multiplicities"),
        "extension_search_bound": (int, "Finite search bound, not a nonexistence certificate"),
    })

    def compute(self, curve, params):
        from sage.all import ZZ, PolynomialRing, GF

        ell, bound = params["ell"], params["extension_bound"]
        if not ZZ(ell).is_prime():
            raise ValueError("ell must be prime")
        if ell == curve.field().characteristic():
            raise AnalysisUnavailable("Characteristic-p torsion needs a finite group-scheme/Dieudonne computation")
        if ell > 65537 or bound > 64:
            raise AnalysisUnavailable("ell exceeds 65537 or extension_bound exceeds 64", "resource_limit")
        polynomial = list(reversed([int(c) for c in l_polynomial(curve, params)]))
        exact, full, least = None, None, None
        if ell == 2 and curve.is_hyperelliptic():
            exact, full, least = hyperelliptic_two_torsion(curve)
            if exact.charpoly("T") != PolynomialRing(GF(2), "T")(polynomial):
                raise ValueError("Branch two-torsion and Frobenius characteristic polynomials disagree")
        rows, semisimple, unipotent, factors = torsion_constraints(polynomial, ell, bound, exact)
        if least is None:
            least = next((row["extension_degree"] for row in rows if row["rank_lower_bound"]), None)
        if full is not None:
            semisimple = full
            while semisimple % ell == 0:
                semisimple //= ell
        elif semisimple is not None and unipotent == 1:
            full = semisimple
        return {"ell": ell, "module_dimension": 2*curve.genus(), "frobenius_factors_mod_ell": factors,
                "rank_information_source": "even branch subsets modulo the full branch set" if exact is not None else "characteristic polynomial; Jordan blocks not assumed",
                "frobenius_matrix": field_matrix(exact) if exact is not None else None,
                "extension_ranks": rows, "base_rational_rank": rows[0]["rank"],
                "base_rank_lower_bound": rows[0]["rank_lower_bound"], "base_rank_upper_bound": rows[0]["rank_upper_bound"],
                "least_nontrivial_extension_degree": least, "semisimple_frobenius_order": semisimple,
                "full_torsion_extension_degree": full,
                "full_torsion_extension_upper_bound": full if full is not None else (semisimple*unipotent if semisimple else None),
                "unipotent_order_upper_bound": unipotent, "extension_search_bound": bound}
