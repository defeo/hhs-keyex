
################################################################################
#
#	A Proof of Concept for key exchange based on isogeny graphs
#
################################################################################

using Nemo
using Primes
using EllipticCurves
using ClassPolynomials
import Base.step #just to use the name

################################################################################
#
#	Things that definitely do not belong here
#
################################################################################


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

function lift(x::FinFieldElem)
	Fl = parent(x)
	@assert Nemo.isprime(order(Fl))
	i = ZZ(0)
	while i != x
		i += 1
	end
	return i
end

################################################################################
#
#	Global parameters
#
################################################################################

p = ZZ(2)^128 + 51
F, _ = FiniteField(p, 1, "u")
FX, X = PolynomialRing(F, "X")

################################################################################
#
#	Walking in the isogeny graph, using modular equations
#
################################################################################


#Compute one step in the ell-isogeny graph, starting from E, in direction v.
#We assume v is given as an integer between 1 and l-1 (Frobenius eigenvalues cannot be zero).
function step(E, ell, v)
	
	#Find neighbors in the ell-isogeny graph
	Phi_ell = ClassicalModularPolynomial(ell, j_invariant(E), X)
	(j1, r1), (j2, r2) = roots(Phi_ell) #'roots' comes from the EllipticCurves module
	
	#There are two simple roots
	@assert r1 == 1
	@assert r2 == 1
	
	#Is j1 correct ?
	phi1 = Isogeny(E, ell, j1)
	K1 = kernel(phi1)
	if is_frobenius_eigenvalue(v, K1, E)
		return image(phi1)
	else
		phi2 = Isogeny(E, ell, j2) #This computation may be removed in the future
		K2 = kernel(phi2)
		@assert is_frobenius_eigenvalue(v, K2, E) #sanity check
		return image(phi2)
	end
end

function next_step(E, ell, v, jprev)
	Phi_ell = ClassicalModularPolynomial(ell, j_invariant(E), X)
	(j1, r1), (j2, r2) = roots(Phi_ell)
	
	#There are two simple roots
	@assert r1 == 1
	@assert r2 == 1
	
	if j1 == jprev
		return j2
	else
		@assert j2 == jprev #sanity check
		return j1
	end
end

function several_steps(E, ell, v, number)
	#we assume number>0
	E1 = step(E, ell, v)
	jprev = j_invariant(E)
	for i = 2:number
		E1, jprev = next_step(E1, ell, v, jprev), j_invariant(E1)
	end
	return E1
end

################################################################################
#
#	Checking Frobenius eigenvalues
#
################################################################################

#Here we check both X- and Y-coordinates to compute the eigenvalue.
#In the case where the Frobenius eigenvalues v1, v2 satisfy v1 != +-v2 mod ell,
#we can just use the X-coordinate and halve the costs.
#Here we assume the eigenvalues are given as integers between 1 and ell-1.

function is_frobenius_eigenvalue(v, K, E)
	_, _, _, a, b = a_invariants(E)
	#Build the coordinate ring
	R = ResidueRing(FX, K(X))
	Xbar = R(X)
	RY, Y = PolynomialRing(R, "Y")
	RR = ResidueRing(RY, Y^2 - Xbar^3 - a * Xbar - b)
	Ybar = RR(Y)
	
	#We assume the elliptic curve is in short Weierstrass form.
	FrobX = power_binexp(Xbar, p)
	FrobY = power_binexp(Xbar^3 + a * Xbar + p, div(p-1, 2)) #we forget the Y factor
	
	Eext = EllipticCurve(RR(a), RR(b))
	Pgeneric = Point(RR(Xbar), Ybar, RR(1), Eext) #the generic point in our kernel
	
	### FAILS ####
	vx, vy, vz = v * Pgeneric
	##############
	
	@assert vz == RR(1) #just in case...
	if (vx == FrobX & vy == Ybar * FrobY)
		return true
	else
		return false
	end
end

################################################################################
#
#	Walking in the isogeny graph, using torsion points
#
################################################################################

function step_torsion(E, ell, Card)
	P = torsionpoint(E, ell, Card)
	return image(Isogeny(E, P))
end


function step_torsion_base_extend(E, ell, Card, r)
	F = base_ring(E)
	p = characteristic(F)
	@assert p == order(F) #here we cannot use non-prime fields
	
	Fext, _ = FiniteField(p, r, "alpha")
	Eext = base_extend(E, Fext)
	Cardext = card_over_extension(Card, p, r)
	
	Pext = torsionpoint(Eext, ell, Cardext)
	
	#Here we can either use the 'Isogeny' function with Pext directly,
	#as above; or we can convert the subgroup back first
	
	Kext = subgroup(Pext, ell)
	K = convert(Kext, F)
	return image(Isogeny(E, K))
end

function several_steps_torsion(E, ell, Card, number)
	#we assume number>=0
	Eprime = E
	for i = 1:number
		Eprime = step_torsion(Eprime, ell, Card)
	end
	return Eprime
end

function several_steps_torsion_base_extend(E, ell, Card, r, number)
	#we assume number>=0
	Eprime = E
	for i = 1:number
		Eprime = step_torsion_base_extend(Eprime, ell, Card, r)
	end
	return Eprime
end

################################################################################
#
#	Key exchange toy example
#
################################################################################


E0 = EllipticCurve(F(1), F(1))
ZZY, Y = PolynomialRing(ZZ, "Y")
frob_poly = Y^2 - ZZ(14234351195508672941)*Y + ZZ(340282366920938463463374607431768211507) #this was computed with Sage
Card = ZZ(340282366920938463449140256236259538567)

#Which primes are we going to use?
#We only want to use primes for which the isogeny graph is a simple cycle, i.e.
#the Frobenius polynomial has exactly two roots mod ell (the frobenius eigenvalues).
#We can use the torsion technique when the two eigenvalues have different multiplicative orders,
#but it is only interesting when the eigenvalue is small enough.
#Of course we have to choose the threshold carefully based on the speed of the algorithms.
#It is likely that this threshold should also depend on ell.

primes_use_step = []
eigenvalues = Dict()
primes_use_step_torsion = []
orders = Dict()

for ell in primes(2, 100) #100 is arbitrary here, but we have classical modular polynomials up to ell=123 only
	Fl, _ = FiniteField(ell, 1, "u")
	FlZ, Z = PolynomialRing(Fl, "Z")
	frob_poly_mod_ell = FlZ(frob_poly)
	rts = roots(frob_poly_mod_ell)
	if length(rts)==2 #this is the only case we are interested in, and excludes the prime 2
		rts = (rts[1][1], rts[2][1])
		order1 = order(rts[1])
		order2 = order(rts[2]) #'order' comes from the EllipticCurves module
		if ((order1 != order2) & (min(order1, order2) < 10)) #5 is aritrary here
			append!(primes_use_step_torsion, [ell])
			orders[ell] = min(order1, order2)
		else
			append!(primes_use_step, ell)
			eigenvalues[ell] = (lift(rts[1]), lift(rts[2]))
		end
	end
end

#Choose bounds for the number of steps to be computed for each prime.
#This determines the size of the key space.
#This choice should be made carefully, as it heavily depends on the cost of the 'step' algorithm for each prime.
#In this toy example, we choose arbitrary bounds.

bounds = Dict()
for ell in primes_use_step
	if ell<20
		bounds[ell] = 5
	else
		bounds[ell] = 1
	end
end
for ell in primes_use_step_torsion
	if orders[ell] == 1
		bounds[ell] = 100
	else
		bounds[ell] = 10
	end
end

#Compute key space size: remember we have no choice for the direction in the torsion case
key_space_size = 1
for ell in primes_use_step
	key_space_size *= (2 * bounds[ell] + 1)
end
for ell in primes_use_step_torsion
	key_space_size *= (bounds[ell] + 1)
end

print("Orders : ")
show(orders)
print("\n")
print("Eigenvalues : ")
show(eigenvalues)
print("\n")
print("Key space size : ")
print(key_space_size)
print("\n")
	
#Key generation
function key_generation(primes_use_step, primes_use_step_torsion, bounds)
	key = Dict()
	for ell in primes_use_step
		key[ell] = rand(-bounds[ell]:bounds[ell])
	end
	for ell in primes_use_step_torsion
		key[ell] = rand(0:bounds[ell])
	end
	return key
end

#Walk in the graph
#In the real world, we should probably compute up to the maximal bound even if the key is smaller
#in order to avoid side-channel attacks, but this is another topic.
function walk(E, key, Card, primes_use_step, primes_use_step_torsion, eigenvalues, orders)
	Eprime = E
	for ell in primes_use_step_torsion
		if key[ell] != 0
			if orders[ell] == 1
				Eprime = several_steps_torsion(Eprime, ell, Card, key[ell])
			else
				Eprime = several_steps_torsion_base_extend(Eprime, ell, Card, orders[ell], key[ell])
			end
		end
	end
	for ell in primes_use_step
		#let's say the first root is the positive direction
		if key[ell] > 0
			Eprime = several_steps(Eprime, ell, eigenvalues[ell][1], key[ell])
		elseif key[ell] < 0
			Eprime = several_steps(Eprime, ell, eigenvalues[ell][2], -key[ell])
		end
	end
	return Eprime
end

#Key exchange proper
function key_exchange_toy_example()
	alice_key = key_generation(primes_use_step, primes_use_step_torsion, bounds)
	bob_key = key_generation(primes_use_step, primes_use_step_torsion, bounds)
	
	alice_curve = walk(E0, alice_key, Card, primes_use_step, primes_use_step_torsion, eigenvalues, orders)
	bob_curve = walk(E0, bob_key, Card, primes_use_step, primes_use_step_torsion, eigenvalues, orders)
	#They publish the curves...
	
	alice_secret = walk(bob_curve, alice_key, Card, primes_use_step, primes_use_step_torsion, eigenvalues, orders)
	bob_secret = walk(alice_curve, bob_key, Card, primes_use_step, primes_use_step_torsion, eigenvalues, orders)
	@assert alice_secret == bob_secret
end

key_exchange_toy_example()



	
