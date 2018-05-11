include("keyex.jl")

p = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769)
Fp, _ = FiniteField(p, 1, "u")

A = Fp(10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171)

E = EllipticCurves.Montgomery(A, Fp(1))

ord = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586670085481138030088461790938201874171652771344144043268298219947026188471598838060)

FpX, _ = PolynomialRing(Fp, "X")

@time params = SystemParams(p, Fp, Dict(
    i => (FiniteField(p, i, "u")[1],
          EllipticCurves.card_over_extension(ord, p, UInt(i)))
    for i in 2:9
), FpX, E, ord, Dict(
    #### 2-way Vélu primes
    3    => VeluPrimeData(409, 1, true, true),
    5    => VeluPrimeData(409, 1, true, true),
    7    => VeluPrimeData(409, 1, true, true),
    11   => VeluPrimeData(409, 1, true, true),
    13   => VeluPrimeData(409, 1, true, true),
    17   => VeluPrimeData(409, 1, true, true),
    103  => VeluPrimeData(409, 1, true, true),
    ##
    19   => VeluPrimeData( 81, 3, true, true),
    ##
    31   => VeluPrimeData( 34, 5, true, true),
    61   => VeluPrimeData( 34, 5, true, true),
    ##
    ##
    29   => VeluPrimeData( 10, 7, true, true),
    71   => VeluPrimeData( 10, 7, true, true),
    ##
    37   => VeluPrimeData(  6, 9, true, true),

    #### 1-way Vélu primes
    523  => VeluPrimeData(409, 1, true, false),
    821  => VeluPrimeData(409, 1, true, false),
    ##
    947  => VeluPrimeData(409, 1, false, true),
    1723 => VeluPrimeData(409, 1, false, true),
    ##
    661  => VeluPrimeData( 81, 3, true, false),
    ##
    1013 => VeluPrimeData(200, 2, false, true),
    1181 => VeluPrimeData(200, 2, false, true),
    ##
    1321 => VeluPrimeData( 34, 5, true, false),
    ##
    547  => VeluPrimeData( 10, 7, true, false),
    ##
    881  => VeluPrimeData(  8, 4, false, true),
    ##
    1693 => VeluPrimeData(  6, 9, true, false),
    ##
    431  => VeluPrimeData( 34, 5, false, true),
    1061 => VeluPrimeData( 34, 5, false, true),

), Dict(
    #### Elkies primes
    23  => ElkiesPrimeData(20, 13, 7, load_Atkin(23, Fp)),
    ##
    41  => ElkiesPrimeData(11, 8, 5, load_Atkin(41, Fp)),
    ##
    43  => ElkiesPrimeData(10, 9, 19, load_Atkin(43, Fp)),
    ##
    47  => ElkiesPrimeData( 9, 17, 11, load_Atkin(47, Fp)),
    ##
    73  => ElkiesPrimeData( 6, 30, 17, load_Atkin(73, Fp)),
    ##
    89  => ElkiesPrimeData( 5, 53, 47, load_Atkin(89, Fp)),
    ##
    107 => ElkiesPrimeData( 4, 64, 5, load_Atkin(107, Fp)),
    109 => ElkiesPrimeData( 4, 35, 28, load_Atkin(109, Fp)),
    113 => ElkiesPrimeData( 4, 94, 6, load_Atkin(113, Fp)),
    ##
    131 => ElkiesPrimeData( 3, 16, 90, load_Atkin(131, Fp)),
    151 => ElkiesPrimeData( 3, 37, 102, load_Atkin(151, Fp)),
    ##
    157 => ElkiesPrimeData( 2, 136, 15, load_Atkin(157, Fp)),
    163 => ElkiesPrimeData( 2, 97, 42, load_Atkin(163, Fp)),
    167 => ElkiesPrimeData( 2, 72, 109, load_Atkin(167, Fp)),
    191 => ElkiesPrimeData( 2, 117, 111, load_Atkin(191, Fp)),
    193 => ElkiesPrimeData( 2, 152, 113, load_Atkin(193, Fp)),
    197 => ElkiesPrimeData( 2, 190, 169, load_Atkin(197, Fp)),
    223 => ElkiesPrimeData( 2, 30, 52, load_Atkin(223, Fp)),
    229 => ElkiesPrimeData( 2, 167, 181, load_Atkin(229, Fp)),
    ##
    241 => ElkiesPrimeData(1, 92, 55, load_Atkin(241, Fp)),
    251 => ElkiesPrimeData(1, 66, 19, load_Atkin(251, Fp)),
    257 => ElkiesPrimeData(1, 211, 95, load_Atkin(257, Fp)),
    277 => ElkiesPrimeData(1, 165, 47, load_Atkin(277, Fp)),
    283 => ElkiesPrimeData(1, 253, 217, load_Atkin(283, Fp)),
    293 => ElkiesPrimeData(1, 219, 99, load_Atkin(293, Fp)),
    307 => ElkiesPrimeData(1, 64, 283, load_Atkin(307, Fp)),
    317 => ElkiesPrimeData(1, 310, 136, load_Atkin(317, Fp)),
    349 => ElkiesPrimeData(1, 79, 53, load_Atkin(349, Fp)),
    359 => ElkiesPrimeData(1, 20, 341, load_Atkin(359, Fp)),
    
    #### Unused Elkies data
    # 3   => ElkiesPrimeData(1, 1, 2, load_Atkin(3, Fp)),
    # 5   => ElkiesPrimeData(1, 1, 4, load_Atkin(5, Fp)),
    # 7   => ElkiesPrimeData(1, 1, 6, load_Atkin(7, Fp)),
    # 11  => ElkiesPrimeData(1, 1, 10, load_Atkin(11, Fp)),
    # 13  => ElkiesPrimeData(1, 1, 12, load_Atkin(13, Fp)),
    # 17  => ElkiesPrimeData(1, 1, 16, load_Atkin(17, Fp)),
    # 19  => ElkiesPrimeData(1, 11, 12, load_Atkin(19, Fp)),
    # 29  => ElkiesPrimeData(1, 25, 22, load_Atkin(29, Fp)),
    # 31  => ElkiesPrimeData(1, 8, 27, load_Atkin(31, Fp)),
    # 37  => ElkiesPrimeData(1, 7, 21, load_Atkin(37, Fp)),
    # 61  => ElkiesPrimeData(1, 20, 3, load_Atkin(61, Fp)),
    # 71  => ElkiesPrimeData(1, 32, 51, load_Atkin(71, Fp)),
    # 103 => ElkiesPrimeData(1, 1, 102, load_Atkin(103, Fp)),
    # ##
    # 373 => ElkiesPrimeData(1, 283, 344, load_Atkin(373, Fp)),
    # 383 => ElkiesPrimeData(1, 46, 27, load_Atkin(383, Fp)),
    # 401 => ElkiesPrimeData(1, 113, 10, load_Atkin(401, Fp)),
    # 421 => ElkiesPrimeData(1, 197, 7, load_Atkin(421, Fp)),
    # 433 => ElkiesPrimeData(1, 351, 34, load_Atkin(433, Fp)),
    # 439 => ElkiesPrimeData(1, 102, 87, load_Atkin(439, Fp)),
    # 443 => ElkiesPrimeData(1, 350, 354, load_Atkin(443, Fp)),
    # 449 => ElkiesPrimeData(1, 399, 329, load_Atkin(449, Fp)),
    # 457 => ElkiesPrimeData(1, 452, 59, load_Atkin(457, Fp)),
    # 467 => ElkiesPrimeData(1, 281, 182, load_Atkin(467, Fp)),
))

## Clean garbage
gc()

## Measure keyspace size
keyspace = BigInt(1)
for (l, d) in params.velu_primes
    keyspace *= 1 + (d.use_left + d.use_right) * d.max_steps
end
for (l, d) in params.elkies_primes
    keyspace *= 1 + 2*d.max_steps
end
println("Keyspace size: ", Int(floor(log2(keyspace))), " bits")

## Some tests on Vélu primes
w = [(Int(l), d.use_right ? right : left, 1) for (l, d) in params.velu_primes]

@assert all([(j_invariant(VeluWalk(E, ell, right, 1, params))
              ≠ j_invariant(VeluWalk(E, ell, left, 1, params)))
             for ell in [19, 31, 61, 29, 71, 37]])

@assert Walk(params, w) == Walk(params, shuffle(w))

## Some tests on Elkies primes
w = [(Int(l), left, 1) for (l, d) in params.elkies_primes]
#@assert Walk(params, w) == Walk(params, shuffle(w))

## The real deal
wv = [(Int(l), d.use_left ? left : right, d.max_steps) for (l, d) in params.velu_primes]
@time Walk(params, wv)
@time begin
    for (l, d) in params.elkies_primes
        Walk(params, [(Int(l), left, d.max_steps)])
    end
end
