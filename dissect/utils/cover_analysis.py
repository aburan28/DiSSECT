"""Verified maps and local ramification on the normalization of curve models."""

from dissect.utils.curve_analysis import function_field, field_value, AnalysisUnavailable


def lcm_degrees(degrees):
    from math import gcd

    result = 1
    for value in degrees:
        result = result*value//gcd(result, value)
    return result


def fiber_record(L, place, different):
    x = L.base_field().gen()
    uniformizer = place.local_uniformizer()
    above = sorted(L.places_above(place), key=str)
    points = []
    for P in above:
        ramification = int(L(uniformizer).valuation(P))
        points.append({
            "place": str(P), "degree": int(P.degree()),
            "ramification_index": ramification,
            "residue_degree": int(P.relative_degree()),
            "different_exponent": int(different.get(P, 0)),
            "wild_different_excess": int(different.get(P, 0))-ramification+1,
            "lower_ramification_break": int(different.get(P, 0))-1 if L.degree() == 2 and ramification == 2 else None,
            "wild": int(ramification % L.characteristic() == 0),
        })
    degree_sum = sum(P["ramification_index"]*P["residue_degree"] for P in points)
    if degree_sum != L.degree():
        raise ValueError("A cover fiber fails the fundamental degree identity")
    return {
        "target_place": str(place), "target_degree": int(place.degree()),
        "at_infinity": int(x.valuation(place) < 0),
        "points": points, "degree_sum": degree_sum,
        "rational_source_points": sum(P["degree"] == 1 for P in points),
        "discriminant_exponent": sum(P["different_exponent"]*P["residue_degree"] for P in points),
    }


def natural_cover(curve, params, L=None):
    L = function_field(curve) if L is None else L
    divisor = L.different()
    different = divisor.dict()
    branch = sorted({P.place_below() for P in different}, key=lambda P: (int(P.degree()), str(P)))
    fibers = [fiber_record(L, Q, different) for Q in branch]
    degree, genus = int(L.degree()), int(L.genus())
    expected = 2*genus-2+2*degree
    if divisor.degree() != expected:
        raise ValueError("The different divisor fails Riemann-Hurwitz")
    rational_base = sorted(L.base_field().places(1), key=str)
    sample = rational_base[:params.get("max_rational_fibers", 64)]
    rational_fibers = [fiber_record(L, Q, different) for Q in sample]
    ramified = [point for fiber in fibers for point in fiber["points"] if point["different_exponent"]]
    monodromy, galois, deck_order, galois_extension = None, None, None, None
    if degree == 2:
        monodromy, galois, deck_order, galois_extension = "C2", 1, 2, 1
    elif degree % int(L.characteristic()) and all(L.polynomial()[i] == 0 for i in range(1, degree)):
        from math import gcd

        q = int(curve.q())
        monodromy, deck_order = "C%d" % degree, gcd(degree, q-1)
        galois = int(deck_order == degree)
        galois_extension = 1
        while pow(q, galois_extension, degree) != 1:
            galois_extension += 1
    return {
        "cover_degree": degree, "target_genus": 0,
        "base_coordinate": str(L.base_field().gen()),
        "extension_coordinate": str(L.gen()),
        "defining_polynomial": str(L.polynomial()), "separable": 1,
        "different_degree": int(divisor.degree()), "riemann_hurwitz_verified": 1,
        "branch_closed_points": len(branch),
        "branch_geometric_points": sum(int(Q.degree()) for Q in branch),
        "ramified_geometric_points": sum(P["degree"] for P in ramified),
        "wild_ramified_geometric_points": sum(P["degree"] for P in ramified if P["wild"]),
        "tame": int(not any(P["wild"] for P in ramified)),
        "branch_fibers": fibers, "rational_fibers": rational_fibers,
        "infinity_fiber": fiber_record(L, L.base_field().place_infinite(), different),
        "branch_splitting_extension_degree": lcm_degrees([int(Q.degree()) for Q in branch]),
        "ramification_splitting_extension_degree": lcm_degrees([P["degree"] for P in ramified]),
        "rational_fibers_complete": int(len(sample) == len(rational_base)),
        "rational_curve_points_from_fibers": sum(F["rational_source_points"] for F in rational_fibers) if len(sample) == len(rational_base) else None,
        "discriminant_degree": sum(F["target_degree"]*F["discriminant_exponent"] for F in fibers),
        "geometric_monodromy_group": monodromy, "galois": galois,
        "deck_group_order": deck_order,
        "galois_constant_extension_degree": galois_extension,
    }


def rational_function(description, L, max_exponent=256):
    """Structured rational expression in [base generator, extension generator]."""
    from dissect.utils.higher_genus_curve import field_element

    field = L.constant_base_field()
    if not isinstance(description, dict) or "raw" in description or "poly" in description:
        return L(field_element(description, field))
    x, y = L(L.base_field().gen()), L.gen()

    def polynomial(terms):
        value = L.zero()
        for term in terms:
            powers = term["powers"]
            if len(powers) != 2 or any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in powers):
                raise ValueError("Map monomials require two nonnegative integer powers")
            if max(powers) > max_exponent:
                raise AnalysisUnavailable("A map exponent exceeds max_map_degree", "resource_limit")
            value += field_element(term["coeff"], field)*x**powers[0]*y**powers[1]
        return value

    numerator = polynomial(description["numerator"])
    denominator = polynomial(description.get("denominator", [{"coeff": 1, "powers": [0, 0]}]))
    if denominator == 0:
        raise ValueError("Map denominator vanishes in the source function field")
    return numerator / denominator


def mapped_cover_ramification(source, target, morphism, u, degree, rational_target=False):
    """Different of a verified separable map, using pullback of differential divisors.

    Ideal generators distinguish different target places above the same base
    place. This works for maps to positive-genus targets as well as P^1.
    """
    from sage.all import PolynomialRing

    base = target if rational_target else target.base_field()
    field = source.constant_base_field()
    pullbacks = {}

    def pullback(Q):
        if Q in pullbacks:
            return pullbacks[Q]
        below = Q if rational_target else Q.place_below()
        image_uniformizer = morphism(target(below.local_uniformizer()))
        candidates = image_uniformizer.divisor_of_zeros().support()
        generators = [morphism(target(g)) for g in Q.prime_ideal().gens()]
        points = []
        for P in candidates:
            ramification = min(int(value.valuation(P)) for value in generators if value)
            if ramification > 0:
                points.append((P, ramification))
        if sum(e*int(P.degree()) for P, e in points) != degree*int(Q.degree()):
            raise ValueError("Pullback of a target place fails the degree identity")
        pullbacks[Q] = points
        return points

    # div(phi^*dt) - phi^*div(dt) is the different of the separable map.
    different = dict(u.differential().divisor().dict())
    for Q, coefficient in target(base.gen()).differential().divisor().dict().items():
        for P, ramification in pullback(Q):
            different[P] = different.get(P, 0)-coefficient*ramification
    different = {P: int(n) for P, n in different.items() if n}
    if any(n < 0 for n in different.values()):
        raise ValueError("A separable map produced a negative different exponent")

    def target_below(P):
        if u.valuation(P) < 0:
            Q0 = base.place_infinite()
        else:
            residue, _, reduce = P.residue_field()
            value = reduce(u)
            orbit, conjugate = [], value
            while conjugate not in orbit:
                orbit.append(conjugate)
                conjugate = conjugate**int(field.cardinality())
            ring = PolynomialRing(residue, "T")
            polynomial = ring.one()
            for conjugate in orbit:
                polynomial *= ring.gen()-conjugate
            constants = {reduce(source(c)): c for c in field}
            coefficients = [constants[c] for c in polynomial]
            prime = sum(c*base.gen()**i for i, c in enumerate(coefficients))
            Q0 = base.maximal_order().ideal(prime).place()
        candidates = [Q0] if rational_target else target.places_above(Q0)
        matches = [Q for Q in candidates if any(place == P for place, _ in pullback(Q))]
        if len(matches) != 1:
            raise ValueError("Could not identify a unique target place below a ramified point")
        return matches[0]

    branch = sorted({target_below(P) for P in different}, key=lambda Q: (int(Q.degree()), str(Q)))
    fibers = []
    for Q in branch:
        points = []
        for P, ramification in pullback(Q):
            residue_degree, remainder = divmod(int(P.degree()), int(Q.degree()))
            if remainder:
                raise ValueError("A residue-field degree is not integral")
            points.append({"place": str(P), "degree": int(P.degree()),
                           "ramification_index": ramification, "residue_degree": residue_degree,
                           "different_exponent": different.get(P, 0),
                           "wild_different_excess": different.get(P, 0)-ramification+1,
                           "lower_ramification_break": different.get(P, 0)-1 if degree == 2 and ramification == 2 else None,
                           "wild": int(ramification % int(field.characteristic()) == 0)})
        fibers.append({"target_place": str(Q), "target_degree": int(Q.degree()), "points": points,
                       "degree_sum": sum(P["ramification_index"]*P["residue_degree"] for P in points)})
    different_degree = sum(n*int(P.degree()) for P, n in different.items())
    target_genus = 0 if rational_target else int(target.genus())
    if different_degree != 2*int(source.genus())-2-degree*(2*target_genus-2):
        raise ValueError("The map different fails Riemann-Hurwitz")
    return {"different_degree": different_degree, "riemann_hurwitz_verified": 1,
            "branch_closed_points": len(branch), "branch_geometric_points": sum(int(Q.degree()) for Q in branch),
            "branch_fibers": fibers,
            "tame": int(not any(P["wild"] for F in fibers for P in F["points"])),
            "wild_ramified_geometric_points": sum(P["degree"] for F in fibers for P in F["points"] if P["wild"])}


def declared_maps(curve, params):
    from sage.all import FunctionField, PolynomialRing, HyperellipticCurve
    from dissect.utils.higher_genus_curve import field_element

    maps = curve.analysis().get("maps", [])
    if len(maps) > params.get("max_maps", 16):
        raise AnalysisUnavailable("Declared maps exceed max_maps", "resource_limit")
    L = function_field(curve)
    output = []
    for specification in maps:
        name = specification.get("name", "unnamed")
        try:
            bound = params.get("max_map_degree", 256)
            u = rational_function(specification["x"], L, bound)
            poles = int(u.divisor_of_poles().degree())
            if not poles:
                raise ValueError("A cover coordinate must be nonconstant")
            target = specification.get("target", {"form": "ProjectiveLine"})
            if target["form"] == "ProjectiveLine":
                target_genus, target_degree = 0, 1
                v = None
                target_field = FunctionField(curve.field(), "t")
                morphism = target_field.hom([u])
            elif target["form"] == "Hyperelliptic":
                ring = PolynomialRing(curve.field(), "t")
                f = ring([field_element(c, curve.field()) for c in target["params"]["f"]])
                h = ring([field_element(c, curve.field()) for c in target["params"].get("h", [0])])
                model = HyperellipticCurve(f, h, check_squarefree=True)
                target_genus, target_degree = int(model.genus()), 2
                v = rational_function(specification["y"], L, bound)
                if v**2 + h(u)*v != f(u):
                    raise ValueError("Images do not satisfy the target equation")
                K = FunctionField(curve.field(), "t")
                ring = PolynomialRing(K, "s")
                target_field = K.extension(ring.gen()**2+h(K.gen())*ring.gen()-f(K.gen()), "s")
                morphism = target_field.hom([v, u])
            elif target["form"] == "SmoothPlane":
                from dissect.utils.higher_genus_curve import HigherGenusCurve

                model = HigherGenusCurve({"name": name+"-target", "category": "map_target",
                                          "form": "SmoothPlane", "field": curve.field_description(),
                                          "params": target["params"]})
                F = model.plane_polynomial()
                v = rational_function(specification["y"], L, bound)
                value = L.zero()
                for powers, coefficient in F.dict().items():
                    value += coefficient*u**powers[0]*v**powers[1]
                if value:
                    raise ValueError("Images do not satisfy F(x,y,1)=0 on the target")
                K = FunctionField(curve.field(), "t")
                ring = PolynomialRing(K, "s")
                equation = ring.zero()
                for powers, coefficient in F.dict().items():
                    equation += coefficient*K.gen()**powers[0]*ring.gen()**powers[1]
                if not equation.derivative():
                    raise AnalysisUnavailable("Target x/z projection is inseparable; choose another plane coordinate system")
                target_genus, target_degree = model.genus(), int(equation.degree())
                target_field = K.extension(equation, "s")
                morphism = target_field.hom([v, u])
            elif target["form"] == "CyclicTrigonal":
                from dissect.utils.higher_genus_curve import HigherGenusCurve

                model = HigherGenusCurve({"name": name+"-target", "category": "map_target",
                                          "form": "CyclicTrigonal", "field": curve.field_description(),
                                          "params": target["params"]})
                v = rational_function(specification["y"], L, bound)
                if v**3 != model.cover_polynomial()(u):
                    raise ValueError("Images do not satisfy the target cubic equation")
                target_field = model.cover_function_field()
                target_genus, target_degree = model.genus(), 3
                morphism = target_field.hom([v, u])
            else:
                raise AnalysisUnavailable("Declared target form is not supported")
            degree, remainder = divmod(poles, target_degree)
            if remainder or degree < 1:
                raise ValueError("Map pole degree is incompatible with the target coordinate")
            if degree > bound:
                raise AnalysisUnavailable("Map degree exceeds max_map_degree", "resource_limit")
            separable = bool(u.derivative())
            inseparable_degree = 1
            if not separable:
                roots = [u] + ([v] if v is not None else [])
                while all(not value.derivative() for value in roots):
                    roots = [value.nth_root(int(curve.field().characteristic())) for value in roots]
                    inseparable_degree *= int(curve.field().characteristic())
            different_degree = 2*curve.genus()-2-degree*(2*target_genus-2) if separable else None
            if different_degree is not None and different_degree < 0:
                raise ValueError("Map violates Riemann-Hurwitz")
            local = mapped_cover_ramification(L, target_field, morphism, u, degree,
                                              rational_target=target["form"] == "ProjectiveLine") if separable else None
            output.append({
                "name": name, "status": "ok", "relation_verified": 1,
                "target_model": target,
                "x_image": str(u), "y_image": str(v) if v is not None else None,
                "degree": degree, "target_genus": target_genus,
                "separable": int(separable), "inseparable_degree": inseparable_degree,
                "separable_degree": degree//inseparable_degree,
                "expected_different_degree": different_degree,
                "local_ramification_computed": int(separable), "ramification": local,
                "reason": None if separable else "The differential different is not assigned to an inseparable cover",
            })
        except (ValueError, KeyError, ZeroDivisionError) as error:
            output.append({"name": name, "status": "invalid_map", "reason": str(error)})
        except AnalysisUnavailable as error:
            output.append({"name": name, "status": error.status, "reason": str(error)})
    for result in output:
        result["source_function_field"] = {"base_coordinate": str(L.base_field().gen()),
                                            "extension_coordinate": str(L.gen()),
                                            "equation": str(L.polynomial())}
    return output


def plane_gonality_cover(curve):
    """Projection through a rational smooth point gives a degree d-1 map."""
    from sage.all import matrix, identity_matrix, FunctionField, PolynomialRing

    points = curve.curve().rational_points()
    if not points:
        return None
    field = curve.field()
    polynomial = curve.plane_polynomial()
    d = int(polynomial.degree())
    basis = identity_matrix(field, 3).columns()
    for point in points:
        for a in basis:
            for b in basis:
                change = matrix(field, [a, list(point), b]).transpose()
                if not change.det():
                    continue
                K = FunctionField(field, "u")
                ring = PolynomialRing(K, "v")
                u, v = K.gen(), ring.gen()
                images = [change[i, 0]*u + change[i, 1]*v + change[i, 2] for i in range(3)]
                equation = ring.zero()
                for powers, coefficient in polynomial.dict().items():
                    term = ring(coefficient)
                    for value, power in zip(images, powers):
                        term *= value**power
                    equation += term
                if equation.degree() != d-1 or not equation.derivative():
                    continue
                L = K.extension(equation, "v")
                if L.genus() != curve.genus():
                    raise ValueError("Projection model changed the geometric genus")
                return {"field": L, "center": [field_value(c) for c in point],
                        "coordinate_change": [[field_value(c) for c in row] for row in change.rows()]}
    raise AnalysisUnavailable("Could not find a separating projection through a rational point")
