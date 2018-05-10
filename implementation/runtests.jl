include("keyex.jl")

# p = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769)
# Fp, _ = FiniteField(p, 1, "u")

# A = Fp(10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171)

# E1 = EllipticCurves.Montgomery(A, Fp(1))
# Card = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586670085481138030088461790938201874171652771344144043268298219947026188471598838060)
# j = j_invariant(E1)

# FpX, X = PolynomialRing(Fp, "X")

p = ZZ(101)
Fp, _ = FiniteField(p, 1, "u")
FpX, X = PolynomialRing(Fp, "X")
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
    7  => ElkiesPrimeData(1, 1, 3, load_Atkin(7, Fp)),
    23 => ElkiesPrimeData(1, 16, 2, load_Atkin(23, Fp)),
))

[
    [ElkiesWalk(E, 3, right, i, params) for i in 1:10],
    [ElkiesWalk(E, 7, left, i, params) for i in 1:10],
    [j_invariant(VeluWalk(M, 3, right, i, params)) for i in 1:10],
    [j_invariant(VeluWalk(M, 7, right, i, params)) for i in 1:10],
    [j_invariant(VeluWalk(M, 61, right, i, params)) for i in 1:10],
    [j_invariant(VeluWalk(M, 67, right, i, params)) for i in 1:10],
    [j_invariant(VeluWalk(M, 409, left, i, params)) for i in 1:10],
]

Walk(params, Dict(
    3  => (left, 5),
    7  => (right, 5),
    23 => (left, 3),
    61 => (right, 2),
    67 => (right, 2),
))
