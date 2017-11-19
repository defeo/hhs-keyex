##############################################################################
#
#		Benchmarks for the HHS key exchange
#
##############################################################################

using Nemo
using BenchmarkTools
using EllipticCurves
using Primes
using ClassPolynomials

suite = BenchmarkGroup()

#We want to measure performance of several key steps of the algorithm:
# * With modular equations:
#	-> Evaluation of the modular polynomial
#	-> Root finding
#	-> Kernel computation
#	-> Computation of the Frobenius endomorphism on the kernel
#	-> Determination of the eigenvalue
#
# * With torsion points:
#	-> Scalar multiplication

suite["modular"] = BenchmarkGroup([])
suite["modular"]["evaluation"] = BenchmarkGroup([])
suite["modular"]["roots"] = BenchmarkGroup([])

suite["direction"] = BenchmarkGroup([])
suite["direction"]["kernel"] = BenchmarkGroup([])
suite["direction"]["eigenvalue"] = BenchmarkGroup([])

suite["torsion"] = BenchmarkGroup([])
suite["torsion"]["weierstrass"] = BenchmarkGroup([])
suite["torsion"]["montgomery"] = BenchmarkGroup([])

##############################################################################
#
#		Parameters
#
##############################################################################

# As a first approximation, we choose a 512-bit prime. This can be readjusted in the future if
# the class group turns out to be too small.
# p is chosen such that p = -1 mod ell for primes ell in [3, 380].
# For maximum accuracy, we will benchmark our algorithms with the actual candidate curve.
# However, we will only benchmark for primes ell with ell < 100, since Atkin modular equations
# can't be used yet.

p = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769)
Fp, _ = FiniteField(p, 1, "u")
a = Fp(11116464848863953015404579872882615582435721470072701808037226672236659358145861235530764279830938006566253824939581765595814774014481950202045501019370057)
b = Fp(7929994832424543515682621940590901949591050994041239823498965175036014605541969152503725727477980546149790133670177545414851226849642510122057033563808087)

E1 = EllipticCurve(a, b)
b, E2 = has_montgomery(E1)
@assert b
#E2 is a curve in Montgomery form, isomorphic to E1 over Fp

Card = ZZ(12037340738208845034383383978222801137092029451270197923071397735408251586670085481138030088461790938201874171652771344144043268298219947026188471598838060)
j = j_invariant(E1)

FpX, X = PolynomialRing(Fp, "X")

##############################################################################
#
#		Benchmarks
#
##############################################################################


for ell in primes(3, 20)
	print("Defining benchmarks for the prime ")
	print(ell)
	print("\n")

	# Evaluation of modular polynomials
	suite["modular"]["evaluation"][(ell, "Classical")] = @benchmarkable ClassicalModularPolynomial($ell, j, X)
	suite["modular"]["evaluation"][(ell, "Atkin")] = @benchmarkable AtkinModularPolynomial($ell, j, X)
	
	# Root finding
	Phi_ell = ClassicalModularPolynomial(ell, j, X)
	suite["modular"]["roots"][ell] = @benchmarkable roots($Phi_ell)
	
	rts = roots(Phi_ell)
	if length(rts) == 2 #make further benchmarks
		(j1, r1), (j2, r2) = rts
		
		#Kernel
		suite["direction"]["kernel"][ell] = @benchmarkable Isogeny(E1, $ell, $j1)
		
		phi1 = Isogeny(E1, ell, j1)
		K1 = kernel(phi1) #this is not a computation
		
		for v in [1, div(ell-1, 2), ell - 1]
			suite["direction"]["eigenvalue"][(ell, v)] = @benchmarkable is_frobenius_eigenvalue($v, $K1, E1)
		end
	end
	
	# Torsion points: it is a good approximation to compute p^d times a random point
	for r = 1:5
		Fext, _ = FiniteField(p, r, "alpha")
		E1ext = base_extend(E1, Fext)
		E2ext = base_extend(E2, Fext)
		suite["torsion"]["weierstrass"][("degree", r)] = @benchmarkable p^($r) * rand($E1ext)    #non-projective law on weierstrass curve
		suite["torsion"]["montgomery"][("degree", r)] = @benchmarkable p^($r) * randXZ($E2ext)  #projective law on montgomery curve
	end
end


################################################################################
#
#	Functions to be benched
#
################################################################################

#Here we check both X- and Y-coordinates to compute the eigenvalue.
#In the case where the Frobenius eigenvalues v1, v2 satisfy v1 != +-v2 mod ell,
#we can just use the X-coordinate and halve the costs.
#Here we assume the eigenvalues are given as integers between 1 and ell-1.

function is_frobenius_eigenvalue(v, K, E)
	_, _, _, a, b = a_invariants(E)
	#Build the coordinate ring
	R = ResidueRing(FpX, K(X))
	Xbar = R(X)
	RY, Y = PolynomialRing(R, "Y")
	RR = ResidueRing(RY, Y^2 - Xbar^3 - a * Xbar - b)
	Ybar = RR(Y)
	
	#We assume the elliptic curve is in short Weierstrass form.
	FrobX = power_binexp(Xbar, p)
	FrobY = power_binexp(Xbar^3 + a * Xbar + b, div(p-1, 2)) #we forget the Y factor
	
	Eext = EllipticCurve(RR(a), RR(b))
	Pgeneric = Point(RR(Xbar), Ybar, RR(1), Eext) #the generic point in our kernel
	
	vx, vy, vz = coordinates(projective_scalar_mul(Pgeneric, v))
	if (vx == FrobX * vz) & (vy == Ybar * FrobY * vz)
		return true
	else
		return false
	end
end

function power_binexp(x, q)
	y = x
	for b in bin(q)[2:end]
		y = y * y
		if (b == '1')
			y = y * x
		end
	end
	return y
end


################################################################################
#
#	Tuning and running
#
################################################################################

print("Tuning benchmarks...\n")

tune!(suite, verbose = true, seconds = 0.1)

print("Running...\n")

results = run(suite, verbose = true, seconds = 0.1)

showall(results)



