from math import gcd
from fractions import Fraction

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import l_polynomial, extension_polynomial, AnalysisUnavailable


def honda_tate_factors(polynomial, prime, field_degree):
    """Simple isogeny factors, including the Schur index at p and real places."""
    from sage.all import NumberField

    blocks = []
    for factor, exponent in polynomial.factor():
        if factor.degree() == 1:
            schur = 2  # Rational real Weil numbers have a real local invariant 1/2.
            local = [{"place": "real", "invariant": [1, 2]}]
        else:
            K = NumberField(factor, "a")
            a = K.gen()
            schur = 2 if K.signature()[0] else 1
            local = [{"place": "real", "invariant": [1, 2]}] if K.signature()[0] else []
            for ideal in K.primes_above(prime):
                valuation = int(a.valuation(ideal))
                residue_degree = int(ideal.residue_class_degree())
                invariant = Fraction(valuation*residue_degree, field_degree) % 1
                schur = schur * invariant.denominator // gcd(schur, invariant.denominator)
                local.append({"place": str(ideal), "valuation": valuation,
                              "residue_degree": residue_degree,
                              "invariant": [invariant.numerator, invariant.denominator]})
        if int(exponent) % schur or int(factor.degree())*schur % 2:
            raise ValueError("Weil polynomial multiplicity fails the Honda-Tate Schur-index condition")
        blocks.append({
            "irreducible_polynomial": [int(c) for c in factor],
            "polynomial_multiplicity": int(exponent), "schur_index": schur,
            "simple_factor_dimension": int(factor.degree())*schur//2,
            "simple_factor_multiplicity": int(exponent)//schur,
            "isotypic_dimension": int(factor.degree())*int(exponent)//2,
            "local_invariants": local,
        })
    return blocks


class JacobianDecompositionTrait(DetailedTrait):
    NAME = "jacobian_decomposition"
    DESCRIPTION = "Honda-Tate isogeny decomposition over the base field and bounded constant extensions."
    INPUT = dict(DetailedTrait.INPUT, max_point_work=(int, "Bound on q^genus"), extension_bound=(int, "Largest splitting-search extension"))
    DEFAULT_PARAMS = dict(DetailedTrait.DEFAULT_PARAMS, max_point_work=[2000000], extension_bound=[3])
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "base_factors": (list, "Simple factors, Schur indices, multiplicities and local invariants"),
        "base_simple": (int, "Jacobian simple up to base-field isogeny"),
        "base_elliptic_factor_dimension": (int, "Total dimension of elliptic isogeny factors"),
        "extension_decompositions": (list, "Exact decompositions in searched extensions"),
        "geometrically_simple": (int, "0 if disproved; null if a bounded search cannot decide"),
        "splitting_witness_extension": (int, "First searched extension proving nonsimplicity"),
        "extension_search_complete": (int, "0: a finite search does not certify absolute simplicity"),
    })

    def compute(self, curve, params):
        from sage.all import PolynomialRing, QQ

        if params["extension_bound"] > 16:
            raise AnalysisUnavailable("extension_bound exceeds 16", "resource_limit")
        ring = PolynomialRing(QQ, "T")
        coefficients = list(reversed([int(c) for c in l_polynomial(curve, params)]))
        p, e = int(curve.field().characteristic()), int(curve.field().degree())
        base = honda_tate_factors(ring(coefficients), p, e)
        simple = len(base) == 1 and base[0]["simple_factor_multiplicity"] == 1
        witness = None if simple else 1
        extensions = []
        for n in range(2, params["extension_bound"]+1):
            factors = honda_tate_factors(ring(extension_polynomial(coefficients, n)), p, e*n)
            extensions.append({"extension_degree": n, "factors": factors})
            if witness is None and (len(factors) > 1 or factors[0]["simple_factor_multiplicity"] > 1):
                witness = n
        return {
            "base_factors": base, "base_simple": int(simple),
            "base_elliptic_factor_dimension": sum(b["simple_factor_multiplicity"] for b in base if b["simple_factor_dimension"] == 1),
            "extension_decompositions": extensions,
            "geometrically_simple": 0 if witness else None,
            "splitting_witness_extension": witness, "extension_search_complete": 0,
        }
