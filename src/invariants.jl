export casimirs

"""
    casimirs(x::Vector, y::Vector) -> (BaseRing, BaseRing)

Compute the Casimir invariants O_n and E_n from cross-ratio coordinates (x_i, y_i).
"""
function casimirs(x::Vector, y::Vector)
    n = length(x)
    if n != length(y)
        error("x and y must have the same length")
    end

    if isodd(n)
        On = prod(x)
        En = prod(y)

        return On, En
    end
   
    if iseven(n)
        On = prod(x[i] for i in 1:n if iseven(i)) + prod(x[i] for i in 1:n if isodd(i))
        En = prod(y[i] for i in 1:n if iseven(i)) + prod(y[i] for i in 1:n if isodd(i))

        return On, En
    end
end

"""
    casimirs(tp::TwistedPolygon) -> (BaseRing, BaseRing)

Compute the Casimir invariants O_n and E_n from the twisted polygon `tp`.
"""
function casimirs(tp::TwistedPolygon)
    xs, ys = cross_ratio_coordinates(tp)
    return casimirs(xs, ys)
end

"""
    monodromy_invariants(tp::TwistedPolygon) -> (BaseRing, BaseRing)

Compute the monodromy invariants Ω_1 and Ω_2 from the twisted polygon `tp`.
"""
function monodromy_invariants(tp::TwistedPolygon)
    M = monodromy(tp)
    Minv = inv(M)
    trM = tr(M)
    detM = det(M)
    xs, ys = cross_ratio_coordinates(tp)
    On, En = casimirs(xs, ys)
    
    Ω_1 = trM^3 / detM
    Ω_2 = tr(Minv)^3 * detM

    return On^2 * En * Ω_1, On * En^2 * Ω_2
end