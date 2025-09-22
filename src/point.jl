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
