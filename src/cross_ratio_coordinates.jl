export cross_ratio_coordinates, pentagram_map, casimirs

"""
    cross_ratio_coordinates(tp::TwistedPolygon)

Compute cross-ratio coordinates (x_i, y_i) for the twisted polygon `tp`.

Definitions used:
  p_i = intersection( (v_i, v_{i+2}), (v_{i+1}, v_{i+3}) )
  x_i = [ v_i, v_{i+1}, p_i, p_{i+1} ]
  y_i = [ v_{i+1}, v_{i+2}, p_i, p_{i+1} ]

Returns two vectors `(xs, ys)` of length n.
"""
function cross_ratio_coordinates(tp::TwistedPolygon)
    verts = tp.vertices
    n = length(verts)
    if n < 5
        error("Need at least 5 vertices to form cross-ratio coordinates.")
    end

    # compute x_i and y_i
    x = Vector{eltype(affine_coords(verts[1]))}(undef, n)
    y = similar(x)
    for i in 1:n
        vi_2 = verts[mod1(i-2, n)]
        vi_1 = verts[mod1(i-1, n)]
        vi   = verts[i]
        vi1 = verts[mod1(i+1, n)]
        vi2  = verts[mod1(i+2, n)]

        # compute x_i
        l1 = line_through(vi_2, vi_1)
        l2 = line_through(vi, vi1)
        p_ip1 = intersection(l1, l2)

        l3 = line_through(vi_2, vi_1)
        l4 = line_through(vi1, vi2)
        p_1p2 = intersection(l3, l4)
    
        x[i] = cross_ratio_on_line(vi_2, vi_1, p_ip1, p_1p2)

        # compute y_i

        l1 = line_through(vi_2,vi_1)
        l2 = line_through(vi1,vi2)
        q_1p1 = intersection(l1, l2)

        l3 = line_through(vi_1, vi)
        l4 = line_through(vi1, vi2)
        q_ip2 = intersection(l3, l4)

        y[i] = cross_ratio_on_line(q_1p1, q_ip2, vi1, vi2)

    end

    return x, y
end

"""
    pentagram_map(x::Vector, y::Vector) -> (Vector, Vector)

Given cross-ratio coordinates, computes the cross-ratio coordinates after one application of the pentagram map.
"""
function pentagram_map(x::Vector, y::Vector)

    # use formula from paper to compute new x and y
    newX = similar(x)
    newY = similar(y)

    for i in eachindex(x)
        newX[i] = x[i] * (1 - x[mod1(i-1, end)] * y[mod1(i-1, end)]) / (1 - x[mod1(i+1, end)] * y[mod1(i+1, end)])
        newY[i] = y[mod1(i+1, end)] * (1 - x[mod1(i+2, end)] * y[mod1(i+2, end)]) / (1 - x[i] * y[i])
    end

    return newX, newY

end

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