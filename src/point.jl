using LinearAlgebra

const BaseRing = Rational
const Point = Vector{BaseRing}

"""
    line_through(a, b) -> Vector{BaseRing}

Return the homogeneous line through points a and b.
"""
line_through(a::AbstractVector{<:BaseRing}, b::AbstractVector{<:BaseRing}) = cross(a, b)

"""
    intersection(l1, l2) -> Vector{BaseRing}

Return the homogeneous intersection of lines l1 and l2.
"""
intersection(l1::AbstractVector{<:BaseRing}, l2::AbstractVector{<:BaseRing}) = cross(l1, l2)


"""
    normalize_point(p::AbstractVector{<:BaseRing}) -> AbstractVector{<:BaseRing}

Normalize a point in homogeneous coordinates so that the last coordinate is 1,
if possible.
"""
function normalize_point(p::AbstractVector{<:BaseRing})
    if p[end] != 0
        return p .// p[end]  # rational normalization if BaseRing <: Integer
    else
        return p
    end
end

"""
    affine_coords(p::AbstractVector{<:BaseRing}) -> Tuple{Float64, Float64}

Convert a projective point [x,y,z] to affine coordinates (x/z, y/z).
"""
function affine_coords(p::AbstractVector{<:BaseRing})
    if p[end] == 0
        error("Point at infinity cannot be projected to affine coordinates")
    end
    return (Float64(p[1] / p[end]), Float64(p[2] / p[end]))
end

"""
    best_coord_index(a_aff::Tuple, b_aff::Tuple) -> Int

Return index (1 or 2) of coordinate with maximal absolute difference between a and b.

This picks a stable coordinate to parametrize the line in the affine chart.
"""
function best_coord_index(a_aff::Tuple, b_aff::Tuple)
    # compare absolute differences in each coordinate
    dx = abs(b_aff[1] - a_aff[1])
    dy = abs(b_aff[2] - a_aff[2])
    return dx >= dy ? 1 : 2
end

function scalar_on_line(a_aff::Tuple, b_aff::Tuple, p_aff::Tuple)
    k = best_coord_index(a_aff, b_aff)
    num = (k == 1) ? (p_aff[1] - a_aff[1]) : (p_aff[2] - a_aff[2])
    den = (k == 1) ? (b_aff[1] - a_aff[1]) : (b_aff[2] - a_aff[2])
    if den == 0
        error("Degenerate line: chosen coordinate difference is zero.")
    end
    return num / den
end

function det3(a::AbstractVector, b::AbstractVector, c::AbstractVector)
    return a[1]*(b[2]*c[3] - b[3]*c[2]) -
           a[2]*(b[1]*c[3] - b[3]*c[1]) +
           a[3]*(b[1]*c[2] - b[2]*c[1])
end

"""
    cross_ratio_on_line(a, b, c, d) -> BaseRing

Note that we use the multiplicative inverse of the "normal" cross ratio formula.
"""
function cross_ratio_on_line(a::AbstractVector, b::AbstractVector,
                             c::AbstractVector, d::AbstractVector; o::AbstractVector = [0,0,1])
    D1 = det3(o,a,b)
    if D1 == 0
        if o == [0,0,1]
            o = [1,0,0]
            return cross_ratio_on_line(a, b, c, d; o=o)
        elseif o == [1,0,0]
            o = [0,1,0]
            return cross_ratio_on_line(a, b, c, d; o=o)
        else
            error("Something is wrong")
        end
    end
    num = det3(o,a,b) * det3(o,c,d)
    den = det3(o,a,c) * det3(o,b,d)

    if den == 0
        error("Choose a different o")
    end
    return num/den
end