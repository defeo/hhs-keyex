###
#
# An implementation of dense bivariate polynomials focusing on
# evaluation, using Julia arrays

using ClassPolynomials
import Nemo: PolynomialRing

struct BiPoly
    coeffs::Matrix{Nemo.fq}
end

function load_Atkin(ell, base)
    _, (F, J) = PolynomialRing(base, ["F", "J"])
    P = AtkinModularPolynomial(ell, F, J)
    degF = maximum(P.exps[1,:])
    degJ = maximum(P.exps[2,:])
    coeffs = Matrix{Nemo.fq}(degF + 1, degJ + 1)
    fill!(coeffs, zero(base))
    for (c, expF, expJ) in zip(P.coeffs, P.exps[1,:], P.exps[2,:])
        coeffs[expF + 1, expJ + 1] = c
    end
    return BiPoly(coeffs)
end

lenJ(P::BiPoly) = size(P.coeffs, 2)
lenF(P::BiPoly) = size(P.coeffs, 1)

function eval_vec(x::T, length::Int) where T
    vec = Vector{T}(length)
    vec[1] = x^0
    for i in 2:length
        vec[i] = vec[i-1] * x
    end
    return vec
end

function diff1(vec::AbstractArray)
    dv = similar(vec)
    dv[1] = zero(vec[1])
    for i in 2:length(vec)
        dv[i] = (i - 1) * vec[i - 1]
    end
    return dv
end

function diff2(vec::AbstractArray)
    dv = similar(vec)
    dv[1] = dv[2] = zero(vec[1])
    for i in 3:length(vec)
        dv[i] = (i - 2) * (i - 1) * vec[i - 2]
    end
    return dv
end

function evalJ(P::BiPoly, Jvec::Vector{T}) where T
    return P.coeffs * Jvec
end

function evalJ(P::BiPoly, J::T) where T
    return evalJ(P, eval_vec(J, lenJ(P)))
end

function evalF(P::BiPoly, Fvec::Vector{T}) where T
    return transpose(Fvec) * P.coeffs
end

function evalF(P::BiPoly, F::T) where T
    return evalF(P, eval_vec(F, lenF(P)))
end

function evalJF(P::BiPoly, Jvec::Vector{T}, Fvec::Vector{U}) where T where U
    return transpose(Fvec) * evalJ(P, Jvec)
end

function evalJF(P::BiPoly, J::T, F::U) where T where U
    return transpose(eval_vec(F, lenF(P))) * evalJ(P, J)
end
