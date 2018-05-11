include("keyex.jl")

######
#
# Small example over GF(101) (h=4)

p = ZZ(101)
Fp, _ = FiniteField(p, 1, "u")
FpX, _ = PolynomialRing(Fp, "X")
E = ShortWeierstrass(Fp(36), Fp(95))
M = Montgomery(Fp(34), Fp(1))
ord = ZZ(120)
params = SystemParams(p, Fp, Dict(
    i => (FiniteField(p, i, "u")[1],
          EllipticCurves.card_over_extension(ord, p, UInt(i)))
    for i in 2:3
), FpX, M, ord, Dict(
    3   => VeluPrimeData(1, 1, true, true),
    7   => VeluPrimeData(1, 1, false, true),
    61  => VeluPrimeData(1, 3, false, true),
    67  => VeluPrimeData(1, 3, false, true),
    409 => VeluPrimeData(1, 3, true, false),
), Dict(
    3  => ElkiesPrimeData(1, 1, 2, load_Atkin(3, Fp)),
    7  => ElkiesPrimeData(1, 4, 6, load_Atkin(7, Fp)),
    23 => ElkiesPrimeData(1, 21, 7, load_Atkin(23, Fp)),
))

@assert ([ElkiesWalk(E, 7, right, i, params) for i in 1:10]
         == [j_invariant(VeluWalk(M, 7, right, i, params)) for i in 1:10])
@assert ElkiesWalk(E, 7, left, 1, params) ≠ ElkiesWalk(E, 7, right, 1, params)
@assert (j_invariant(VeluWalk(M, 3, right, 1, params))
         ≠ j_invariant(VeluWalk(M, 3, left, 1, params)))

@assert Walk(params, [
    (3, left, 5),
    (7, right, 5),
    (23, left, 3),
    (61, right, 2),
    (67, right, 2),
]) == Walk(params, [
    (7, right, 5),
    (3, left, 5),
    (61, right, 2),
    (23, left, 3),
    (67, right, 2),
]) == 38

######
#
# A larger example (128 bits)

p = ZZ(2)^128 + 51
Fp, _ = FiniteField(p, 1, "u")
FpX, _ = PolynomialRing(F, "X")
E = EllipticCurve(Fp(1), Fp(9))
ord = ZZ(340282366920938463477190161499806233971)

# primes_use_step = []
# eigenvalues = Dict()
# primes_use_step_torsion = []
# orders = Dict()

# for ell in primes(2, 100)
#     Fl, _ = FiniteField(ell, 1, "u")
#     FlZ, Z = PolynomialRing(Fl, "Z")
#     frob_poly_mod_ell = FlZ(frob_poly)
#     rts = roots(frob_poly_mod_ell)
#     if length(rts)==2 #this is the only case we are interested in, and excludes the prime 2
# 	rts = (rts[1][1], rts[2][1])
# 	order1 = order(rts[1])
# 	order2 = order(rts[2]) #'order' comes from the EllipticCurves module
# 	if ((order1 != order2) & (min(order1, order2) < 10)) #5 is aritrary here
# 	    append!(primes_use_step_torsion, [ell])
# 	    orders[ell] = min(order1, order2)
# 	else
# 	    append!(primes_use_step, ell)
# 	    eigenvalues[ell] = (lift(rts[1]), lift(rts[2]))
# 	end
#     end
# end
