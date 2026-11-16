"""Exact zeta-series arithmetic and prime-to-p torsion information."""

from math import comb

from dissect.utils.curve_analysis import power_sums, extension_polynomial


def mobius(n):
    sign, prime = 1, 2
    while prime*prime <= n:
        if n % prime == 0:
            n //= prime
            sign = -sign
            if n % prime == 0:
                return 0
        prime += 1
    return -sign if n > 1 else sign


def divisor_counts(coefficients, q, genus, bound):
    """Closed points and rational effective divisors, with Euler-product checks."""
    traces = power_sums(list(reversed(coefficients)), bound)
    points = [None]+[q**n+1-traces[n] for n in range(1, bound+1)]
    closed, effective = [0], [1]
    for n in range(1, bound+1):
        numerator = sum(mobius(d)*points[n//d] for d in range(1, n+1) if n % d == 0)
        count, remainder = divmod(numerator, n)
        if remainder or count < 0:
            raise ValueError("Zeta data did not give nonnegative integral closed-point counts")
        closed.append(count)
        effective.append((q+1)*effective[-1]-(q*effective[-2] if n > 1 else 0)
                         +(coefficients[n] if n < len(coefficients) else 0))
    euler = [1]+[0]*bound
    for degree in range(1, bound+1):
        if not closed[degree]:
            continue
        factor = [comb(closed[degree]+k-1, k) for k in range(bound//degree+1)]
        euler = [sum(euler[n-k*degree]*factor[k] for k in range(n//degree+1)) for n in range(bound+1)]
    if euler != effective:
        raise ValueError("The effective-divisor series failed its closed-point Euler product")
    order = sum(coefficients)
    rr = []
    for n in range(max(1, 2*genus-1), bound+1):
        expected = order*(q**(n+1-genus)-1)//(q-1)
        if effective[n] != expected:
            raise ValueError("The effective-divisor count failed Riemann-Roch")
        rr.append(n)
    if bound >= genus and effective[genus] != order+q*effective[genus-2]:
        raise ValueError("The degree-genus Picard-fiber identity failed")
    return points, closed, effective, rr


def hyperelliptic_two_torsion(curve):
    """Frobenius on even branch subsets modulo the full branch set."""
    from sage.all import GF, matrix, identity_matrix, lcm

    f, h = curve.hyperelliptic_polynomials()
    branch = h*h+4*f
    degrees = [int(P.degree()) for P, _ in branch.factor()]
    if branch.degree() % 2:
        degrees.append(1)
    size = sum(degrees)
    if size != 2*curve.genus()+2:
        raise ValueError("The hyperelliptic branch divisor has the wrong degree")
    permutation, start = [], 0
    for degree in degrees:
        permutation.extend(start+(i+1) % degree for i in range(degree))
        start += degree
    F = matrix(GF(2), size-2)
    for i in range(size-2):
        image = [0]*size
        image[permutation[i]] ^= 1
        image[permutation[-1]] ^= 1
        for j in range(size-2):
            F[j,i] = image[j] ^ image[size-2]
    order = lcm(degrees)
    identity = identity_matrix(GF(2), size-2)
    for prime, _ in order.factor():
        while order % prime == 0 and F**(order//prime) == identity:
            order //= prime
    least = min(int(n) for n in order.divisors() if (F**n-identity).rank() < size-2)
    return F, int(order), least


def polynomial_power_mod(base, exponent, modulus):
    # Sage 9.8's GF(2) PowerMod path aborts on unreduced degree-one bases.
    # Reducing every product also keeps this independent of the NTL backend.
    result = base.parent().one()
    base %= modulus
    while exponent:
        if exponent % 2:
            result = (result*base) % modulus
        base = (base*base) % modulus
        exponent //= 2
    return result


def torsion_constraints(polynomial, ell, bound, exact_matrix=None):
    """Bounds over every Jordan form with the given characteristic polynomial."""
    from sage.all import GF, PolynomialRing, identity_matrix

    ring = PolynomialRing(GF(ell), "T")
    T = ring.gen()
    reduced = ring(polynomial)
    factors = list(reduced.factor())
    rows, semisimple_order = [], None
    for n in range(1, bound+1):
        ell_part, coprime_part = 1, n
        while coprime_part % ell == 0:
            ell_part *= ell
            coprime_part //= ell
        lower = upper = 0
        for factor, multiplicity in factors:
            if polynomial_power_mod(T, n, factor) == 1:
                lower += int(factor.degree())*min(int(multiplicity), ell_part)
                upper += int(factor.degree())*int(multiplicity)
        if upper == reduced.degree() and semisimple_order is None:
            semisimple_order = n
        if exact_matrix is not None:
            rank = int(exact_matrix.nrows()-(exact_matrix**n-identity_matrix(exact_matrix.base_ring(), exact_matrix.nrows())).rank())
            if not lower <= rank <= upper:
                raise ValueError("The branch-subset module contradicts the Frobenius polynomial")
            lower = upper = rank
        extension = extension_polynomial(polynomial, n)
        order = sum(extension)
        valuation, residual = 0, order
        while residual % ell == 0:
            valuation += 1
            residual //= ell
        if valuation < lower or (valuation > 0) != (lower > 0):
            raise ValueError("Torsion ranks contradict the Jacobian extension order")
        upper = min(upper, valuation)
        rows.append({"extension_degree": n, "rank_lower_bound": lower, "rank_upper_bound": upper,
                     "rank": lower if lower == upper else None, "jacobian_order": order,
                     "ell_adic_order_valuation": valuation})
    largest = max(int(m) for _, m in factors)
    unipotent_bound = 1
    while unipotent_bound < largest:
        unipotent_bound *= ell
    return rows, semisimple_order, unipotent_bound, [
        {"polynomial": [int(c) for c in P], "degree": int(P.degree()), "multiplicity": int(m)} for P, m in factors]
