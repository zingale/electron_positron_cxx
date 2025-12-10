import mpmath as mp


def legendre_pair(n, x):
    """Compute P_n(x) and P_{n-1}(x) (Legendre polynomials) using the
    three-term recurrence, in mpmath arithmetic.

    """
    x = mp.mpf(x)
    if n == 0:
        return mp.mpf('1'), mp.mpf('0')  # P_0, P_-1 (dummy)
    if n == 1:
        return x, mp.mpf('1')           # P_1, P_0

    Pnm1 = mp.mpf('1')   # P_0
    Pn   = x             # P_1
    for k in range(1, n):
        # P_{k+1}(x) = ((2k+1)x P_k(x) - k P_{k-1}(x)) / (k+1)
        Pkp1 = ((2*k + 1) * x * Pn - k * Pnm1) / (k + 1)
        Pnm1, Pn = Pn, Pkp1
    # on exit: Pn = P_n, Pnm1 = P_{n-1}
    return Pn, Pnm1


def gauss_legendre_mpmath(n, digits=100):
    """n-point Gauss–Legendre nodes and weights on [-1, 1] with
    ~`digits` decimal digits of precision using mpmath.

    Returns (nodes, weights) as lists of mp.mpf, sorted ascending.

    """
    # working precision: add a guard factor
    mp.mp.dps = digits + 20

    nodes   = [mp.mpf('0')] * n
    weights = [mp.mpf('0')] * n

    half = n // 2
    tol  = mp.mpf(10) ** (-(digits + 5))
    maxit = 50

    for k in range(1, half + 1):
        # initial guess (classical asymptotic formula)
        theta = mp.pi * (4*k - 1) / (4*n + 2)
        x = mp.cos(theta)

        # Newton iteration to refine root
        for _ in range(maxit):
            Pn, Pnm1 = legendre_pair(n, x)
            # (1 - x^2) P_n'(x) = n (P_{n-1}(x) - x P_n(x))
            dPn = n * (Pnm1 - x * Pn) / (1 - x**2)
            dx = Pn / dPn
            x -= dx
            if abs(dx) < tol:
                break

        # NOW recompute Pn, Pn−1 at *final* x
        Pn, Pnm1 = legendre_pair(n, x)
        dPn = n*(Pnm1 - x*Pn)/(1 - x**2)

        # symmetric nodes
        i_left  = k - 1
        i_right = n - k
        nodes[i_left]  = -x
        nodes[i_right] =  x

        # weight: w = 2 / ((1 - x^2) * P_n'(x)^2)
        w = 2 / ((1 - x**2) * dPn**2)
        weights[i_left]  = w
        weights[i_right] = w

    # odd n: central node at x = 0
    if n % 2 == 1:
        mid = n // 2
        x0 = mp.mpf('0')
        Pn, Pnm1 = legendre_pair(n, x0)
        dPn = n * (Pnm1 - x0 * Pn) / (1 - x0**2)
        w0 = 2 / ((1 - x0**2) * dPn**2)
        nodes[mid]   = x0
        weights[mid] = w0

    # nodes are already in ascending order by construction
    return nodes, weights


def gauss_laguerre_mpmath(n, digits=100):
    """n-point Gauss–Laguerre quadrature for ∫₀^∞ f(x)  dx.

    Returns (x, w) where:
      - x[i]  are the nodes  (x > 0)
      - w[i]  are the weights * exp(x)
    computed with ~`digits` decimal digits of precision using mpmath.

    This uses the Golub–Welsch algorithm for the Laguerre polynomials
    L_n(x), which are orthonormal on [0, ∞) with weight e^{-x}.

    """
    # working precision with a guard
    mp.mp.dps = digits + 20

    # Jacobi matrix J for Laguerre polynomials L_k(x) (α = 0):
    #   x L_k(x) = (2k+1) L_k(x) - k L_{k-1}(x) - (k+1)L_{k+1}(x)
    # On the basis {L_0, ..., L_{n-1}} the multiplication-by-x operator
    # is represented by the symmetric tridiagonal matrix:
    #
    #   J[i,i]   = 2*i + 1
    #   J[i,i-1] = J[i-1,i] = -i   (for i >= 1)
    #
    J = mp.matrix(n)
    for i in range(n):
        J[i, i] = 2*i + 1
        if i > 0:
            J[i, i-1] = -mp.mpf(i)
            J[i-1, i] = -mp.mpf(i)

    # symmetric eigendecomposition: eigenvalues -> nodes, eigenvectors -> weights
    vals, vecs = mp.eigsy(J)

    # sort by node value
    pairs = sorted(
        ((vals[i], vecs[:, i]) for i in range(n)),
        key=lambda p: p[0]
    )

    xs = [p[0] for p in pairs]
    # For orthonormal polynomials, Laguerre in this normalization,
    # the Gauss weights are w_i = (first component of eigenvector)^2
    ws = [(p[1][0])**2 for p in pairs]

    # Reduce precision back down for user-facing values
    # and add exp(+x) to the weights
    mp.mp.dps = digits
    xs = [+x for x in xs]  # unary + rounds to current mp.mp.dps
    ws = [+w * mp.exp(x) for x, w in zip(xs, ws)]

    return xs, ws


if __name__ == "__main__":
    n = 800
    digits = 100

    xs, ws = gauss_legendre_mpmath(n, digits=digits)

    # check sum of weights ~ 2
    print("sum(weights) =", mp.nsum(lambda i: ws[int(i)], [0, n-1]))
    print()

    mp.mp.dps = digits

    # print a few points/weights with 100 digits
    for i in range(n):
        print(f"x[{i}] = {xs[i]}")
        print(f"w[{i}] = {ws[i]}\n")

    xs, ws = gauss_laguerre_mpmath(n, digits=digits)

    mp.mp.dps = digits

    print("Laguerre: ")

    # print a few points/weights with 100 digits
    for i in range(n):
        print(f"x[{i}] = {xs[i]}")
        print(f"w[{i}] = {ws[i]}\n")
