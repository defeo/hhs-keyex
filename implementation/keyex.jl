using Nemo
#using Primes
import EllipticCurves: ShortWeierstrass, Montgomery, Point, XZPoint, unsafe_kernelpoly, j_invariant, from_j_invariant, coordinates, projective_scalar_mul, base_extend, xdouble, xadd, xinfinity, isinfinity, normalized

include("bipolys.jl")
include("c_libs.jl")

###

@enum Direction::Int8 left=-1 right=1

Degree = UInt64

struct VeluPrimeData
    max_steps::UInt16
    ev_order::UInt16
    use_left::Bool
    use_right::Bool
end

struct ElkiesPrimeData
    max_steps::UInt16
    ev_left::Degree
    ev_right::Degree
    mod_pol::BiPoly
end

struct SystemParams
    p::Nemo.fmpz
    Fp::Nemo.FqFiniteField
    FpExts::Dict{UInt8, Tuple{Nemo.FqFiniteField, Nemo.fmpz}}
    FpX::Nemo.FqPolyRing
    E0::Montgomery
    order::Nemo.fmpz
    velu_primes::Dict{Degree, VeluPrimeData}
    elkies_primes::Dict{Degree, ElkiesPrimeData}
end

###

"""

Compute `steps` isogeny steps from `E` along the direction
`eigenvalue`, using Elkies' algorithm.

"""
function ElkiesWalk(E::ShortWeierstrass,
                    ell::Integer,
                    dir::Direction,
                    steps::Integer,
                    params::SystemParams)
    j = j_invariant(E)
    if steps > 0
        prime_data = params.elkies_primes[ell]
        Ψ = prime_data.mod_pol
        eigenvalue = dir == right ? prime_data.ev_right : prime_data.ev_left

        # First step
        j_vec = eval_vec(j, lenJ(Ψ))

        # We get two roots
        Ψ_F_j = evalJ(Ψ, j_vec)          # Ψ(F, j)
        roots = my_roots(params.FpX(Ψ_F_j))
        @assert length(roots) == 2
        (f1, r1), (f2, r2) = roots
        @assert r1 == r2 == 1

        # We try our chance with the first root
        f_vec = eval_vec(f1, lenF(Ψ))
        Ψ_f_J = evalF(Ψ, f_vec)          # Ψ(f1, J)
        f_quo, f_rem = divrem(params.FpX(transpose(Ψ_f_J)), gen(FpX) - j)
        @assert f_rem == 0
        roots = my_roots(f_quo)

        # try each root until one works
        j1 = zero(params.Fp)
        φ = one(params.FpX)
        for (j1, r) in roots
            # Elkies' formulas
            j1_vec  = eval_vec(j1, lenJ(Ψ))
            dj_vec  = diff1(j_vec)
            dj1_vec = diff1(j1_vec)
            df_vec  = diff1(f_vec)
            Ψ_df_j  = transpose(df_vec) * Ψ_F_j
            Ψ_df_j1 = evalJF(Ψ, j1_vec, df_vec)
            Ψ_f_dj, Ψ_f_dj1 = Ψ_f_J * hcat(dj_vec, dj1_vec)

            dF  = f1 * Ψ_df_j
            dJ  = j  * Ψ_f_dj
            dF1 = f1 * Ψ_df_j1
            dJ1 = j1 * Ψ_f_dj1 * ell
            
            jj = j1 // (j1 - 1728)
            tmp = ell^2 * dF1 * dJ * E.b // (dJ1 * dF * E.a)
            # Image curve
            a1 = -27 * tmp^2 * jj // 4
            b1 = a1 * tmp
            E1 = ShortWeierstrass(a1, b1)

            # We get the kernel polynomial and check the eigenvalue
            issq, φ = issquare(unsafe_kernelpoly(E, E1, ell)[2])
            if issq
                break
            end
        end
        @assert φ ≠ 1
        # @assert EllipticCurves.divisionpolynomial(E, ell) % φ == 0

        # TODO: handle multiple Atkin roots
        if is_frobenius_eigenvalue(eigenvalue, φ, E, params)
            f = f1
            j = j1
        else
            # If we chose the wrong direction, we get the opposite j-invariant
            f = f2
            Ψ_f_J = evalF(Ψ, f)          # Ψ(f, J)
            f_quo, f_rem = divrem(params.FpX(transpose(Ψ_f_J)), gen(FpX) - j)
            @assert f_rem == 0
            (j, r), = my_roots(f_quo)
            @assert r == 1
        end

        # Finally, we walk the remaining steps
        for step in 2:steps
            Ψ_F_j = evalJ(Ψ, j)
            f_quo, f_rem = divrem(params.FpX(Ψ_F_j), gen(FpX) - f)
            @assert f_rem == 0
            (f, r), = my_roots(params.FpX(f_quo))
            @assert r == 1
            ###
            Ψ_f_J = evalF(Ψ, f)
            f_quo, f_rem = divrem(params.FpX(transpose(Ψ_f_J)), gen(FpX) - j)
            @assert f_rem == 0
            (j, r), = my_roots(f_quo)
            @assert r == 1
        end
    end        
    return j
end


function is_frobenius_eigenvalue(eigenvalue::Integer,
                                 kernel_poly::Nemo.fq_poly,
                                 E::ShortWeierstrass,
                                 params::SystemParams)
    # TODO: use X coordinate only
    R = ResidueRing(params.FpX, kernel_poly(gen(FpX)))
    Xbar = R(gen(FpX))
    RY, Y = PolynomialRing(R, "Y")
    RR = ResidueRing(RY, Y^2 - Xbar^3 - E.a * Xbar - E.b)
    Ybar = RR(Y)

    FrobX = Xbar^BigInt(params.p)
    FrobY = (Xbar^3 + E.a * Xbar + E.b)^BigInt(div(params.p - 1, 2)) #we forget the Y factor
    
    Eext = ShortWeierstrass(RR(E.a), RR(E.b))
    Pgeneric = Point(RR(Xbar), Ybar, RR(1), Eext) #the generic point in our kernel
    vx, vy, vz = coordinates(projective_scalar_mul(Pgeneric, eigenvalue))
    return (vx == FrobX * vz) && (vy == Ybar * FrobY * vz)
end

###

"""

Compute `steps` isogeny steps from `E` along the direction
`eigenvalue`, using Vélu's formulas.

"""
function VeluWalk(E::Montgomery,
                  ell::Integer,
                  dir::Direction,
                  steps::Integer,
                  params::SystemParams)
    if steps > 0
        prime_data = params.velu_primes[ell]
        ext = prime_data.ev_order
        if ext == 1
            Fext, ord = params.Fp, params.order
        else
            Fext, ord = params.FpExts[ext]
        end
        
        # correct if we use points on the twist
        if dir == left
            @assert prime_data.use_left
            legendre = 1
        else
            @assert prime_data.use_right
            legendre = -1
            ord = 2*p^Int(ext) + 2 - ord
        end

        # Start walking
        for step in 1:steps
	    Eext = base_extend(E, Fext)
            cofactor = div(ord, ell)
            P = xinfinity(Eext)
            i = 1
            while isinfinity(P) || !isinfinity(ell*P)
                P = XZPoint(gen(Fext) + i, one(Fext), Eext)
                i += 1
                if (P.X*(P.X^2 + Eext.A*P.X + 1)*Eext.B)^div(order(Fext) - 1, 2) != legendre
                    continue
                end
                P = cofactor*P
            end
            Q, R = P, xdouble(P)
            σ = σ1 = zero(Fext)
            π = one(Fext)
            for i in 1:div(ell-1, 2)
                σ  += Q.X // Q.Z
                σ1 += Q.Z // Q.X
                π  *= Q.X // Q.Z
                Q, R = R, xadd(R, P, Q)
            end
            σ = Nemo.convert(σ, params.Fp)
            σ1 = Nemo.convert(σ1, params.Fp)
            π = Nemo.convert(π, params.Fp)
	    E = Montgomery((6*σ1 - 6*σ + E.A)*π^2, E.B)
        end
    end
    return E
end

###

function Walk(params::SystemParams, key::Array{T}) where T <: Tuple{Integer, Direction, Integer}
    elkies = []
    E = params.E0
    for (ell, dir, steps) in key
        if haskey(params.velu_primes, ell)
            E = VeluWalk(E, ell, dir, steps, params)
        else
            push!(elkies, (ell, dir, steps))
        end
    end
    j = j_invariant(E)
    for (ell, dir, steps) in elkies
        E = from_j_invariant(j)
        j = ElkiesWalk(E, ell, dir, steps, params)
    end
    return j
end
