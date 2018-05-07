
# coding: utf-8

# In[3]:

from sage.arith.all import is_square
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.matrix.berlekamp_massey import berlekamp_massey
from sage.rings.infinity import Infinity


# In[4]:

# The names by which the algorithms of this module may be called
bmss_names = ("BMSS", "bmss")
stark_names = ("Stark", "stark", "Starks", "starks")
algorithm_names = bmss_names + stark_names


# In[5]:

def isogeny_kernel(E1, E2, degree, A, P, S, algorithm="BMSS"):
    r"""
    Compute the kernel polynomial of the normalized rational isogeny between
    ``E1`` and ``E2`` of degree ``degree`` .

    Assuming a rational normalized isogeny of degree ``degree`` exists between
    ``E1`` and
    ``E2``, this function returns the squarefree polynomial vanishing on the
    abscissae of its kernel. If no such isogeny exists, the outcome depends on
    the algorithm:

        - If the BMSS algorithm is used, its output is tested for correctness:
          if the test succeeds an exception is returned, if it fails a random
          polynomial is returned. The test has heuristic failure probability
          exponentially small in ``degree``, see :py:func:`isogeny_BMSS`
          for details.

        - If Stark's algorithm is used, an exception is returned. (This is
          most probably not true. Luca De Feo, Jul 29, 2011)

    Better checks are implemented by
    :py:class:`sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
    and by
    :py:meth:`sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny`,
    which are the preferred ways to create isogenies at the sage command line.

    :param E1: an elliptic curve in short Weierstrass form.
    :param E2: an elliptic curve in short Weierstrass form.
    :param degree: the degree of the isogeny from ``E1`` to ``E2``.
    :param A: an univariate polynomial ring over the base field.
    :param P: a power series ring over the base field.
    :param S: a Laurent series ring over the base field.
    :param algorithm: ``BMSS`` (default) or ``Stark``.

    :returns: the squarefree polynomial vanishing on the abscissae of the kernel
        of the isogeny.
    :rtype: Polynomial

    :raises ZeroDivisionError: when the characteristic is too small to compute
        isogenies of degree ``degree``. See below.
    :raises ValueError: when no isogeny is found.
    :raises ValueError: when the two curves are not in short Weierstrass form.

    .. WARNING::

        If the characteristic of the base field is `p>0`, the algorithms of this
        module only work when a certain bound on ``degree`` is satisfied:

            - for BMSS, ``p > 4 degree - 1``,
            - for Stark's, ``p > 4 degree + 6``.

        If the bound is not satisfied, a ``ZeroDivisionError`` is raised.

    EXAMPLES:

    We can compute isogenies over the rationals::

        sage: from sage.schemes.elliptic_curves.ell_isogeny_char_zero import isogeny_kernel

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: R.<x> = QQ[]
        sage: f = x
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain(); E2
        Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field
        sage: isogeny_kernel(E, E2, 2)
        x

    However beware that BMSS may return a wrong answer::

        sage: E2 = EllipticCurve([0,0,0,16,0])
        sage: isogeny_kernel(E, E2, 2)
        x
        sage: E.isogeny(None, codomain=E2, degree=2)
        Traceback (most recent call last):
        ...
        ValueError: Codomain parameter must be isomorphic to computed codomain isogeny

    While Stark's algorithm is somewhat safer (but see the documentation of
    :py:func:`isogeny_BMSS` for another way of avoiding wrong answers)::

        sage: isogeny_kernel(E, E2, 2, "Stark")
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a rational normalized isogeny of degree 2

    We can also compute isogenies over finite fields::

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: isogeny_kernel(E, E2, 5)
        x^2 + 7*x + 13
        sage: f
        x^2 + 7*x + 13

    and number fields, too::

        sage: R.<x> = QQ[]
        sage: K.<i> = NumberField(x^2 + 1)
        sage: E = EllipticCurve(K, [0,0,0,1,0])
        sage: E2 = EllipticCurve(K, [0,0,0,16,0])
        sage: isogeny_kernel(E, E2, 4, algorithm="stark")
        x^3 + x

    """
    p = E1.base_ring().characteristic()
    X = A.gen()

    # BMSS algorithm
    if algorithm in bmss_names or algorithm is None:
        if p > 0 and p <= 4*degree-1:
            raise ZeroDivisionError("BMSS algorithm only works for "
                                    "characteristic 0 or greater than 4*degree - 1.")

        # BMSS returns the denominator of the x-component of a normalized
        # isogeny of degree at most degree,
        # or a non-sense polynomial if no such isogeny exists.
        ker_poly = isogeny_BMSS(E1, E2, degree, P)
        # Here we check that the isogeny has EXACTLY the required degree.
        if ker_poly.degree() != degree-1:
            raise ValueError("The two curves are not linked by a rational "
                             "normalized isogeny of degree %s" % degree)
        # Here we check that it probably is an isogeny and not any random polynomial.
        # See the discussion in the documentation.
        even_part = ker_poly(X).gcd(E1.two_division_polynomial()(X))
        odd_part, _ = ker_poly(X).quo_rem(even_part)
        check, sqodd_part = is_square(odd_part, root=True)
        if not check:
            raise ValueError("The two curves are not linked by a rational "
                             "normalized isogeny of degree %s" % degree)
        return even_part * sqodd_part

    # Stark's algorithm
    elif algorithm in stark_names:
        if p > 0 and p <= 4*degree+6:
            raise ZeroDivisionError("Stark's algorithm only works for characteristic 0 or greater than 4*degree + 6.")

        # Stark's algorithm returns the denominator of the x-component of the
        # isogeny, we get its squarefree part
        ker_poly = isogeny_Stark(E1, E2, degree, A, S)
        even_part = ker_poly(X).gcd(E1.two_division_polynomial()(X))
        odd_part, _ = ker_poly(X).quo_rem(even_part)
        check, sqodd_part = is_square(odd_part, root=True)
        if not check:
            raise RuntimeError("Stark's algorithm returned an unexpected result. Please report this bug.")
        return even_part * sqodd_part

    else:
        raise ValueError("Unknown algorithm '%s'." % algorithm)


# In[6]:

def isogeny_Stark(E1, E2, degree, A, S):
    r"""
    Compute the kernel of the degree ``degree`` isogeny between ``E1`` and ``E2`` via
    Stark's algorithm.  There must be a degree ``degree``, rational, separable,
    normalized isogeny from ``E1`` to ``E2``.

    :param E1: an elliptic curve in short Weierstrass form.
    :param E2: an elliptic curve in short Weierstrass form.
    :param degree: the degree of the isogeny from ``E1`` to ``E2``.
    :param A: an univariate polynomial ring over the base field
    :param S: a Laurent series ring over the base field

    :returns: the polynomial vanishing on the abscissae of the kernel
        of the isogeny. Notice that the abscissa of points of order greater than
        2 is contributed twice in the polynomial (once from the point and once
        from its opposite).
    :rtype: Polynomial

    :raises ZeroDivisionError: when the characteristic is smaller or equal to
        ``4*degree + 6`` and not 0.
    :raises ValueError: when no isogeny is found.
    :raises ValueError: when the two curves are not in short Weierstrass form.

    ALGORITHM:

    This function uses Stark's Algorithm as presented in section 6.2 of
    [BMSS]_.

    .. note::

       As published there, the algorithm is incorrect, and a correct
       version (with slightly different notation) can be found in
       [M09]_.  The algorithm originates in [S72]_

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_isogeny_char_zero import isogeny_Stark

        sage: E1 = EllipticCurve(GF(97), [52, 31])
        sage: R.<x> = GF(97)[]; f = x^5 + 67*x^4 + 13*x^3 + 35*x^2 + 77*x + 69
        sage: phi = EllipticCurveIsogeny(E1, f)
        sage: E2 = phi.codomain()
        sage: isogeny_Stark(E1, E2, 11)
        x^10 + 37*x^9 + 53*x^8 + 66*x^7 + 66*x^6 + 17*x^5 + 57*x^4 + 6*x^3 + 89*x^2 + 53*x + 8

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: isogeny_Stark(E, E2, 5)
        x^4 + 14*x^3 + x^2 + 34*x + 21
        sage: f**2
        x^4 + 14*x^3 + x^2 + 34*x + 21

    """

    K = E1.base_field()
    x = A.gen()

    try:
        wp1 = E1.weierstrass_p(prec=4*degree+4)  #BMSS claim 2*degree is enough, but it is not M09
        wp2 = E2.weierstrass_p(prec=4*degree+4)
    except NotImplementedError:
        raise ZeroDivisionError("Stark's algorithm only works for characteristic 0 or greater than 4*degree + 6.")

    # viewed them as power series in Z = z^2
    Z = S.gen()
    pe1 = 1/Z
    pe2 = 1/Z
    for i in xrange(2*degree+1):
        pe1 += wp1[2*i] * Z**i
        pe2 += wp2[2*i] * Z**i
    pe1 = pe1.add_bigoh(2*degree+2)
    pe2 = pe2.add_bigoh(2*degree+2)

    #print 'wps = ',pe1
    #print 'wps2 = ',pe2

    n = 1
    q = [R(1), R(0)]
    #p = [R(0), R(1)]
    T = pe2

    while ( q[n].degree() < (degree-1) ):
        #print 'n=', n

        n += 1
        a_n = 0
        r = -T.valuation()
        while (0 <= r):
            t_r = T[-r]
            #print '    r=',r
            #print '    t_r=',t_r
            #print '    T=',T
            a_n = a_n + t_r * x**r
            T = T - t_r*pe1**r
            r = -T.valuation()


        q_n = a_n*q[n-1] + q[n-2]
        q.append(q_n)
        #p_n = a_n*p[n-1] + q[n-2]
        #p.append(p_n)

        if n == degree + 1 or T == 0:
            if T == 0 or T.valuation() < 2:
                raise ValueError("The two curves are not linked by a rational normalized isogeny of degree %s" % degree)
            #print 'breaks here'
            break

        T = 1/T
        #print '  a_n=', a_n
        #print '  q_n=', q_n
        #print '  p_n=', p_n
        #print '  T = ', T

    qn = q[n]
    #pn= p[n]
    #print 'final  T = ', T
    #print '  f =', pn/qn

    qn = (1/qn.leading_coefficient())*qn
    #pn = (1/qn.leading_coefficient())*pn

    return qn


# In[7]:

def isogeny_BMSS(E1, E2, degree, P):
    r"""
    Compute the kernel of the normalized isogeny between ``E1`` and ``E2`` via
    the BMSS algorithm.  There must be a rational, separable,
    normalized isogeny of degree at most ``degree`` from ``E1`` to ``E2``.

    :param E1: an elliptic curve in the Weierstrass form `y^2 = f(x)`.
    :param E2: an elliptic curve in the Weierstrass form `y^2 = f(x)`.
    :param degree: a bound on the degree of the isogeny from E1 to E2.
    :param P: a power series ring over the base field.

    :returns: the polynomial vanishing on the abscissae of the kernel
        of the isogeny, if the isogeny exists, or a random looking polynomial
        otherwise. Notice that if the isogeny exists, the abscissa of points of
        order greater than
        2 is contributed twice in the polynomial (once from the point and once
        from its opposite).
    :rtype: Polynomial

    :raises ZeroDivisionError: when the characteristic is smaller or equal to
        ``4*degree - 1`` and not 0.
    :raises ValueError: when no isogeny is found.
    :raises ValueError: when the two curves are not in the Weierstrass form
        `y^2 = f(x)`.

    .. WARNING::

        If no rational normalized isogeny of degree at most ``degree``
        exists between ``E1`` and ``E2``, this function silently fails by
        returning a somewhat random polynomial of degree ``degree - 1``.
        A check on this result must be performed in order to tell it apart
        from the successfull case.

        One easy test is to check that the result is the product of a square
        polynomial and of a factor of the 2-division polynomial. Obviously the
        probability that a random polynomial satisfies this condition depends
        upon what one means by "random", not to mention that the result of this
        algorithm is obviously not random. However it is fair to say that the
        failure probability of this test decreases exponentially with
        ``degree``.

        Other, more constraining, test are implemented by
        :py:class:`sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`,
        and by
        :py:meth:`sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny`,
        which are the preferred ways to create isogenies at the sage command
        line.

    .. note::

        The BMSS algorithm differs from Stark's algorithm in that it doesn't
        need the degree of the isogeny as input. However, for efficiency (and to
        insure termination in characteristic `0`), this implementation takes a
        bound on the degree of the isogeny as input.

        The key property is that, given two elliptic
        curves defined over a field of characteristic `p`, there is at most one
        normalized isogeny of degree less than `(p+1)/4` (or of any degree in
        case `p=0`), and the BMSS algorithm will eventually find it.

        To prove the uniqueness, suppose that `\phi_1` and `\phi_2` are both
        normalized, then by [Sil86]_, III.5.2

        .. math::

            (\phi_1 - \phi_2)^\ast \omega_2 = 0

        so `\phi_1-\phi_2` is not separable ([Sil86]_, II.4.2(c)). This rules
        out the case `p=0`.

        Now suppose that `p>0`. Denote by `d(\phi)` the degree of an isogeny
        `\phi`, and suppose w.l.o.g. that `d(\phi_1)>d(\phi_2)`.
        Because `\phi_1-\phi_2` is unseparable, `d(\phi_1-\phi_2)`
        is at least `p` ([Sil86]_, II.2.12). The degree is a positive definite
        quadratic form on `\mathrm{Hom}(E_1,E_2)` ([Sil86]_, III.6.3), so by
        the triangle inequality ([Sil86]_, V.1.2)

        .. math::

            p \le
            d(\phi_1 - \phi_2) \le
            \left(\sqrt{d(\phi_1)} + \sqrt{d(\phi_2)}\right)^2 \le
            4d(\phi_1) <
            p+1.

        Hence `p=4d(\phi_1)`, contraddicting the fact that `p` is prime. The
        credit for this proof goes to Marco Streng.

        The proof that the algorithm will indeed find the isogeny is in [BMSS]_.
        The only addendum is that the degree of the isogeny only plays
        a role in determining the truncation order of a certain power series
        (see below),
        and that computing beyond this truncation order does not change the
        result.

    ALGORITHM:

    This function uses an unpublished variant of the algorithm presented in
    section 4.3 of [BMSS]_. Let `N/D` be the rational fraction expressing the
    action of the isogeny on the abscissae, we compute the power series
    `T(x) = D(1/x)/N(1/x)`, solution to the equation

    .. math::

       (a_6 x^3 + a_4 x^2 + a_2 x + 1) {T'}^2 =
       (T/x) (\tilde{a}_6 T^3 + \tilde{a}_4 T^2 + \tilde{a}_2 T + 1)

    with initial conditions `T = x + O(x^2)`; then we reconstruct the rational
    fraction using the :py:func:`Berlekamp-Massey algorithm
    <sage.matrix.berlekamp_massey.berlekamp_massey>`.

    EXAMPLES:

    ::

        sage: from sage.schemes.elliptic_curves.ell_isogeny_char_zero import isogeny_BMSS, isogeny_kernel

    Here we don't know the degree of the isogeny, so we give an upper bound of
    8; the algorithm correctly finds the degree 5 normalized isogeny::

        sage: E = EllipticCurve(GF(37), [0,0,0,1,8])
        sage: R.<x> = GF(37)[]
        sage: f = (x + 14) * (x + 30)
        sage: phi = EllipticCurveIsogeny(E, f)
        sage: E2 = phi.codomain()
        sage: kernel = isogeny_BMSS(E, E2, 8)
        sage: kernel
        x^4 + 14*x^3 + x^2 + 34*x + 21
        sage: phi2 = EllipticCurveIsogeny(E, kernel.sqrt())
        sage: phi2.degree()
        5
        sage: phi == phi2
        True

    However, the output may be wrong if no isogeny of degree at most ``degree``
    exists::

        sage: E1 = EllipticCurve(GF(97), [52, 31])
        sage: R.<x> = GF(97)[]; f = x^5 + 67*x^4 + 13*x^3 + 35*x^2 + 77*x + 69
        sage: phi = EllipticCurveIsogeny(E1, f)
        sage: E2 = phi.codomain()
        sage: kernel = isogeny_BMSS(E1, E2, 9)
        sage: kernel
        x^8 + 36*x^7 + 95*x^6 + 84*x^4 + 16*x^3 + 70*x^2 + 12*x + 78
        sage: kernel.is_square()
        False
        sage: isogeny_kernel(E1, E2, 9, algorithm="bmss")
        Traceback (most recent call last):
        ...
        ValueError: The two curves are not linked by a rational normalized isogeny of degree 9

    Also recall the example from :py:func:`isogeny_kernel`, showing that when
    the bound on the degree is very low, the chances a wrong answer passing
    the test in :py:func:`isogeny_kernel` are high::

        sage: E = EllipticCurve(QQ, [0,0,0,1,0])
        sage: E2 = EllipticCurve([0,0,0,16,0])
        sage: f = isogeny_BMSS(E, E2, 2); f
        x
        sage: f.divides(E.division_polynomial(2))
        True
        sage: isogeny_BMSS(E, E2, 3)
        x^2 + 4/3

    Test for :trac:`11095`::

        sage: from sage.schemes.elliptic_curves.ell_isogeny_char_zero import *
        sage: E = EllipticCurve([-1,0])
        sage: E2 = EllipticCurve([-3^4,0])
        sage: E.isogeny(kernel=None, codomain=E2, degree=9)
        Isogeny of degree 9 from Elliptic Curve defined by y^2 = x^3 - x over Rational Field to Elliptic Curve defined by y^2 = x^3 - 81*x over Rational Field

    """

    (a1, a2, a3, a4, a6) = E1.a_invariants()
    (b1, b2, b3, b4, b6) = E2.a_invariants()

    if a1 != 0 or a3 != 0 or b1 != 0 or b3 != 0:
        raise ValueError("Curves must have a model of the form y^2 = f(x).")

    K = E1.base_field()
    x = P.gen()

    G = a6*x**3 + a4*x**2 + a2*x + 1
    H = b6*x**3 + b4*x**2 + b2*x + 1

    # solve the differential equation
    # G(x) T'^2 = (T/x) H(T)
    T = _BMSS_diffeq(G, H, 2*degree + 1)

    # We recover the rational fraction using the relation
    # T == x * D.reverse() / N.reverse()
    U = T.shift(-1)
    N = berlekamp_massey(U.padded_list())
    D = (U * P(N.reverse())).truncate().reverse()

    # If the points of abscissa 0 are in the kernel,
    # correct the degree of D
    gap = N.degree() - D.degree() - 1
    if (gap > 0): D = D.shift(gap)

    return D


# In[8]:

def _BMSS_diffeq(G, H, prec=None):
    r"""
    Compute a power-series solution to the differential equation used by
    the BMSS algorithm. The differential equation is

    .. math::  P(T,x) = (T/x) H(T) - G(x) {T'}^2 = 0

    with initial conditions `T = O(x)`. The output is truncated to
    precision ``prec``.

    :param G: a power series with non-zero constant coefficient.
    :param H: a power series with non-zero constant coefficient.
    :param prec: (optional) an integer denoting the truncation order of the solution.
        If not given, it defaults to the common precision of G and H, or to the
        default precision of the parent ring.

    :returns: the solution to the differential equation.
    :rtype: Power series

    :raises ZeroDivisionError: when the characteristic is smaller than
        ``2*prec-2`` and not 0.
    :raises ValueError: when ``G`` or ``H`` is not invertible.

    ALGORITHM:

    It is a Newton iteration. After some substitution we get

    .. math::

       T_0 = ax + O(x^2)\\
       T_{i+1} = T_i + T_i' \sqrt{G} \sqrt{x} \int \frac{k_i(x)}{2\sqrt{x}}
       + O(x^{2^i+1}),

    where

    .. math::

       k_i(x) = \frac{P(T_i, x)} {{T_i'}^2 G \sqrt{G}}.

    EXAMPLES::

       sage: from sage.schemes.elliptic_curves.ell_isogeny_char_zero import _BMSS_diffeq

       sage: P.<x> = PowerSeriesRing(QQ, 'x', 20)
       sage: G = H = 0
       sage: while not G.is_unit(): G = P.random_element() ^ 2
       sage: while not H.is_unit(): H = P.random_element()
       sage: T = _BMSS_diffeq(G, H)
       sage: H(T) * T.shift(-1) - G * T.derivative()^2
       O(x^19)

    """

    # The power series ring
    R = G.parent()
    p = R.base_ring().characteristic()

    if not G.is_unit() or not H.is_unit():
        raise ValueError("No unique solution: arguments must be "
                         "invertible power series.")

    if prec is None:
        prec = G.common_prec(H)

    if prec >= Infinity:
        prec = R.default_prec()

    if 0 < p and p <= 2*prec-3:
        raise ZeroDivisionError("Characteristic must be greater than 2*prec - 3 = " +
                                str(2*prec-3) +
                                " in order to compute a solution to precision 'prec'.")

    # the precision to which the solution T is known
    d = 1
    # 1/(T'^2 G) to precision d
    diffT2G = G.O(1)/(H.O(1)**2)
    # Sqrt(G) and Sqrt(1/G) to precision d
    sqG = G.O(1).sqrt(); invsqG = ~sqG

    T = (H.O(1) / G.O(1)).shift(1)

    while d < prec-1:
        # update diffT, sqG and invsqG to precision d
        # (nothing changes in the first iteration)
        diffT2G = diffT2G * (2 - G * T.derivative()**2 * diffT2G)
        sqG = (sqG + G * invsqG * (2 - sqG * invsqG)) / 2
        invsqG = invsqG * (2 - sqG * invsqG)

        # double the current precision
        d = min(2*d, prec-1)
        diffT2G = R(list(diffT2G)).O(d)
        sqG = R(list(sqG)).O(d)
        invsqG = R(list(invsqG)).O(d)
        T = R(list(T)).O(d+1)

        k = (H(T) * T.shift(-1) - G * T.derivative()**2) * diffT2G * invsqG
        # K = 1/2 sqrt(x) integral(k/sqrt(x))
        K =  R([c/(2*i+1) for (i, c) in enumerate(k)]).O(k.prec()).shift(1)

        # update the solution
        T += T.derivative() * sqG * K

    return T


# In[ ]:

def compute_p_1(E1, E2, l, Phi):
    r"""
    Compute the sum of the roots of the kernel polynomial of an isogeny of degree ''l'' linking ''E1'' to ''E2''.
    
    :param E1: an Elliptic curve in Weierstrass form y^2 = x^3 + A*x + B
    :param E2: an Elliptic curve in Weierstrass form y^2 = x^3 + A*x + B over the same field
    :param l: the isogeny degree
    :param Phi: the modular polynomial of degree ''l''
    
    :return: the sum of the roots of the polynomial defining the kernel of this isogeny, as given by Schoof's formulae []
    :rtype: Field Element
    
    There must be a rational, separable, cyclic isogeny of degree ''l'' between ''E1'' and ''E2''
    """
    
    #Field, invariants and polynomials
    K=E1.base_field()
    _, _, _, A1, B1 = E1.a_invariants()
    _, _, _, A2, B2 = E2.a_invariants()
    j1 = E1.j_invariant()
    j2 = E2.j_invariant()
    X, Y = Phi.parent().gens()
    Phi_X = Phi.derivative(X)
    Phi_Y = Phi.derivative(Y)
    Phi_XX = Phi_X.derivative(X)
    Phi_XY = Phi_X.derivative(Y)
    Phi_YY = Phi_Y.derivative(Y)
    
    #Building Eisenstein series E_4, E_6
    E4_1 = - K(48)*A1
    E6_1 = - K(864)*B1
    
    E4_2 = - K(48)*A2
    E6_2 = - K(864)*B2
    
    a_1 = E6_1 / E4_1
    b_1 = E4_1**2 / E6_1
    
    a_2 = E6_2 / E4_2
    b_2 = E4_2**2 / E6_2
     
    #Schoof's formulae
    #Computing the power series dj/dq
    dj1 = - a_1 * j1
    dj2 = - a_2 * j2

    #Computing d2j1/dj1-l*d2j2/dj2
    numerator = - (dj1**2)*Phi_XX(j1,j2) + 2*l*dj1*dj2*Phi_XY(j1,j2) + (l**2)*(dj2**2)*Phi_YY(j1,j2)
    denominator = j2*Phi_X(j1,j2)
    C_1 = numerator / denominator 
    
    #Computing E2_1-l*E2_2
    
    C_2 = 6 * (C_1 + K(1/2)*b_1 + K(2/3)*a_1 - K(l/2)*b_2 - K(2*l/3)*a_2)
    
    #Computing p1
    p_1 = K(1/12)*l*C_2
    
    return p_1


# In[ ]:



