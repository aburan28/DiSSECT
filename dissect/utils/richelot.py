"""The fifteen Richelot kernels and their odd-characteristic quotient targets."""

from math import prod

from dissect.utils.curve_analysis import field_value, l_polynomial, extension_polynomial
from dissect.utils.cover_analysis import lcm_degrees
from dissect.utils.curve_automorphisms import extension_field
from dissect.utils.projective_configurations import encoded, value_key, field_description, frobenius_configuration


def pair_partitions(indices):
    if not indices:
        yield ()
        return
    first = indices[0]
    for i in range(1, len(indices)):
        for rest in pair_partitions(indices[1:i]+indices[i+1:]):
            yield ((first, indices[i]),)+rest


def partition_action(partitions, permutation):
    canonical = lambda partition: tuple(sorted(tuple(sorted(pair)) for pair in partition))
    lookup = {canonical(partition): i for i, partition in enumerate(partitions)}
    action = [lookup[canonical([(permutation[a], permutation[b]) for a,b in partition])] for partition in partitions]
    seen, cycles = set(), []
    for i in range(len(action)):
        if i in seen:
            continue
        cycle, j = [], i
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = action[j]
        cycles.append(cycle)
    degrees = [0]*len(action)
    for cycle in cycles:
        for i in cycle:
            degrees[i] = len(cycle)
    return {"permutation": action, "cycles": cycles, "orbit_degrees": sorted(map(len, cycles)),
            "rational_over_base": sum(len(cycle) == 1 for cycle in cycles)}, degrees


def genus_two_l_polynomial(polynomial):
    """Independent direct counts for the small Richelot verification fields."""
    from sage.all import PolynomialRing

    field, counts = polynomial.base_ring(), []
    q = int(field.cardinality())
    for degree in (1, 2):
        extension, embedding = extension_field(field, degree)
        ring = PolynomialRing(extension, "x")
        F = ring([embedding(c) for c in polynomial])
        roots = lambda value: 1 if not value else (2 if value.is_square() else 0)
        counts.append(sum(roots(F(x)) for x in extension)+roots(F[6]))
    trace, trace_two = q+1-counts[0], q*q+1-counts[1]
    middle, remainder = divmod(trace*trace-trace_two, 2)
    if remainder:
        raise ValueError("Richelot target counts gave a nonintegral zeta coefficient")
    return [1, -trace, middle, -q*trace, q*q]


def richelot_targets(curve, params):
    from sage.all import PolynomialRing, matrix

    f, h = curve.hyperelliptic_polynomials()
    source = f+h*h/4
    degrees = [int(P.degree()) for P,_ in source.factor()]
    if source.degree() % 2:
        degrees.append(1)
    if sum(degrees) != 6:
        raise ValueError("Richelot construction requires six distinct branch points")
    splitting = lcm_degrees(degrees)
    partitions = list(pair_partitions(tuple(range(6))))
    result = {"geometric_kernel_count": 15, "branch_splitting_extension_degree": splitting,
              "kernel_order": 4, "polarized_isogeny_degree": 4, "kernel_maximally_isotropic": 1,
              "targets_complete": 0}
    if int(curve.q())**splitting > params["max_field_size"]:
        permutation, start = [], 0
        for degree in degrees:
            permutation.extend(start+(i+1) % degree for i in range(degree))
            start += degree
        action, definition_degrees = partition_action(partitions, permutation)
        return dict(result, kernel_frobenius=action, base_field_kernel_count=action["rational_over_base"],
                    kernel_partitions=[{"pairs": [list(pair) for pair in P], "kernel_definition_degree": definition_degrees[i]} for i,P in enumerate(partitions)],
                    status="partial", reason="Branch splitting field exceeds max_field_size; abstract kernels and their Frobenius action remain exact")
    field, embedding = extension_field(curve.field(), splitting)
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    source = ring([embedding(c) for c in source])
    points = [(root, field.one()) for root in source.roots(multiplicities=False)]
    if source.degree() % 2:
        points.append((field.one(), field.zero()))
    points.sort(key=lambda P: tuple(value_key(c) for c in P))
    if len(points) != 6:
        raise ValueError("The declared branch field did not split the six points")
    branch_action = frobenius_configuration(points, int(curve.q()))
    action, definition_degrees = partition_action(partitions, branch_action["permutation"])
    verify = int(field.cardinality())**2 <= params["max_target_point_work"]
    expected = list(reversed(extension_polynomial(list(reversed([int(c) for c in l_polynomial(curve, params)])), splitting))) if verify else None
    double_ring = PolynomialRing(field, names=("u", "v"))
    u, v = double_ring.gens()
    targets, records, split_count, verified = [], [], 0, 0
    for i, partition in enumerate(partitions):
        factors = [(points[a][1]*x-points[a][0])*(points[b][1]*x-points[b][0]) for a,b in partition]
        factors[0] *= source.leading_coefficient()/prod(factors).leading_coefficient()
        if prod(factors) != source:
            raise ValueError("The Richelot quadratic factors do not multiply to the source equation")
        determinant = matrix(field, [[P[j] for j in range(3)] for P in factors]).det()
        brackets = [factors[(j+1)%3].derivative()*factors[(j+2)%3]-factors[(j+1)%3]*factors[(j+2)%3].derivative() for j in range(3)]
        if sum(P(u)*B(v) for P,B in zip(factors, brackets))+determinant*(u-v)**2:
            raise ValueError("The Richelot fundamental bilinear identity failed")
        record = {"pairs": [list(pair) for pair in partition], "kernel_definition_degree": definition_degrees[i],
                  "quadratic_factors": [[field_value(c) for c in P] for P in factors],
                  "determinant": field_value(determinant), "brackets": [[field_value(c) for c in P] for P in brackets],
                  "fundamental_identity_verified": 1}
        if determinant:
            target = prod(brackets)/determinant
            if target.degree() not in (5,6) or not target.is_squarefree():
                raise ValueError("A nondegenerate Richelot target must be a smooth genus-two equation")
            target_l = genus_two_l_polynomial(target) if verify else None
            if verify and target_l != expected:
                raise ValueError("Richelot target and source Jacobians have different Weil polynomials")
            verified += int(verify)
            record.update(target_kind="genus_two_jacobian", target_f=[field_value(c) for c in target], target_h=[],
                          target_l_polynomial=target_l, source_target_weil_polynomial_verified=1 if verify else None)
            targets.append(record["target_f"])
        else:
            split_count += 1
            record.update(target_kind="product_of_elliptic_curves", target_components_computed=0)
        records.append(record)
    result.update(branch_points=[encoded(P) for P in points], branch_frobenius=branch_action,
                  target_field=field_description(field, embedding, curve.field()), kernel_frobenius=action,
                  base_field_kernel_count=action["rational_over_base"], kernel_partitions=records,
                  genus_two_target_count=len(targets), split_product_kernel_count=split_count,
                  distinct_target_equation_count=len({str(P) for P in targets}), targets_complete=1,
                  point_count_verified_targets=verified,
                  target_descent_scope="Targets are explicit over the branch splitting field; a smaller kernel field does not certify descent of this chosen equation")
    return result
