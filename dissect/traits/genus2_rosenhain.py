from itertools import permutations

from dissect.traits.detailed_trait import DetailedTrait
from dissect.utils.curve_analysis import check_budget, field_value, AnalysisUnavailable
from dissect.utils.cover_analysis import lcm_degrees
from dissect.utils.curve_automorphisms import extension_field
from dissect.utils.projective_configurations import encoded, value_key, frobenius_configuration, field_description


class Genus2RosenhainTrait(DetailedTrait):
    NAME = "genus2_rosenhain"
    DESCRIPTION = "Six-point branch configurations, Rosenhain orbits and rational Jacobian two-torsion."
    GENUS = 2
    SUPPORTED_CURVE_TYPES = ("hyperelliptic",)
    OUTPUT = dict(DetailedTrait.OUTPUT, **{
        "branch_orbit_degrees": (list, "Closed branch-point degrees over the input field"),
        "branch_splitting_extension_degree": (int, "Least constant extension splitting all six branch points"),
        "rosenhain_model_over_base": (int, "All six branch points are base-field rational; a scalar twist is allowed"),
        "rational_two_torsion_rank": (int, "Dimension of J[2](F_q), using the even-subset quotient representation"),
        "splitting_field": (dict, "Field basis and embedding used for configurations"),
        "branch_points": (list, "Ordered homogeneous P1 branch coordinates"),
        "branch_frobenius": (dict, "Explicit Frobenius permutation and cycles"),
        "rosenhain_configurations": (list, "Distinct unordered lambda triples after ordered 0,1,infinity normalizations"),
        "configuration_count": (int, "Number of distinct unordered triples"),
        "ordered_normalizations_tested": (int, "Number of distinct ordered triples chosen from six points"),
        "canonical_configuration": (list, "Smallest triple in the reported field basis; does not distinguish twists"),
        "configuration_complete": (int, "Whether the complete six-point normalization orbit was computed"),
    })

    def compute(self, curve, params):
        from sage.all import PolynomialRing, matrix, GF

        check_budget(curve, params)
        if curve.field().characteristic() == 2:
            raise AnalysisUnavailable("Classical six-point Rosenhain configurations require odd characteristic")
        f, h = curve.hyperelliptic_polynomials()
        branch = h**2+4*f
        degrees = [int(factor.degree()) for factor, _ in branch.factor()]
        if branch.degree() % 2:
            degrees.append(1)
        if sum(degrees) != 6:
            raise ValueError("A smooth odd-characteristic genus-two double cover must have six branch points")
        splitting_degree = lcm_degrees(degrees)
        permutation = []
        start = 0
        for degree in degrees:
            permutation.extend(start+(i+1) % degree for i in range(degree))
            start += degree
        # Fixed vectors in {sum(v)=0}/<ones> can satisfy (F-1)v=ones,
        # so merely counting irreducible factors is insufficient in general.
        equations = matrix(GF(2), 7, 7)
        for i, image in enumerate(permutation):
            equations[image, i] += 1
            equations[i, i] += 1
            equations[i, 6] = 1
            equations[6, i] = 1
        rank = int(equations.right_kernel().dimension())-1
        result = {"branch_orbit_degrees": sorted(degrees),
                  "branch_splitting_extension_degree": splitting_degree,
                  "rosenhain_model_over_base": int(splitting_degree == 1),
                  "rational_two_torsion_rank": rank, "configuration_complete": 0}
        if int(curve.q())**splitting_degree > params["max_field_size"]:
            return dict(result, status="partial", reason="Branch splitting field exceeds max_field_size; orbit degrees and two-torsion rank remain exact")
        field, embedding = extension_field(curve.field(), splitting_degree)
        polynomial = PolynomialRing(field, "x")([embedding(c) for c in branch])
        points = [(root, field.one()) for root in polynomial.roots(multiplicities=False)]
        if branch.degree() % 2:
            points.append((field.one(), field.zero()))
        if len(points) != 6:
            raise ValueError("The branch splitting field did not contain six distinct points")
        points.sort(key=lambda P: tuple(value_key(c) for c in P))
        configurations = {}
        for i, j, k in permutations(range(6), 3):
            a, b, c = points[i], points[j], points[k]
            determinant = lambda P, Q: P[1]*Q[0]-P[0]*Q[1]
            values = []
            for index, P in enumerate(points):
                if index not in (i, j, k):
                    values.append(determinant(a, P)*determinant(c, b)/(determinant(c, P)*determinant(a, b)))
            values.sort(key=value_key)
            configurations[tuple(value_key(v) for v in values)] = values
        keys = sorted(configurations)
        result.update(splitting_field=field_description(field, embedding, curve.field()),
                      branch_points=[encoded(P) for P in points],
                      branch_frobenius=frobenius_configuration(points, int(curve.q())),
                      rosenhain_configurations=[[field_value(v) for v in configurations[key]] for key in keys],
                      configuration_count=len(keys), ordered_normalizations_tested=120,
                      canonical_configuration=[field_value(v) for v in configurations[keys[0]]],
                      configuration_complete=1)
        return result
