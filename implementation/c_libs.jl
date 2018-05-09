import Nemo

"""
    Flint's own root-finding algorithm, as called via *_poly_factor_equal_deg
    """
function flint_roots(x::Nemo.fq_poly)
    R = Nemo.parent(x)
    F = Nemo.base_ring(R)
    fac = Nemo.fq_poly_factor(F)
    ccall((:fq_poly_factor_equal_deg, :libflint), Void, 
          (Ref{Nemo.fq_poly_factor}, Ref{Nemo.fq_poly}, Int,
           Ref{Nemo.FqFiniteField}), fac, x, 1, F)
    res = Dict{Nemo.fq, Int}()
    for i in 1:fac.num
        f = R()
        ccall((:fq_poly_factor_get_poly, :libflint), Void,
              (Ref{Nemo.fq_poly}, Ref{Nemo.fq_poly_factor}, Int,
               Ref{Nemo.FqFiniteField}), f, fac, i-1, F)
        d = unsafe_load(fac.exp,i)
        res[-Nemo.coeff(f,0)] = d
    end
    return res
end

function my_roots(x::Nemo.fq_poly)
    R = Nemo.parent(x)
    X = Nemo.gen(R)
    F = Nemo.base_ring(R)
    q = Nemo.order(F)
    Xq = Nemo.powmod(X, q, x)
    lin = Nemo.gcd(x, Xq - X)
    if Nemo.degree(lin) > 0
        return flint_roots(lin)
    else
        return Dict()
    end
end
