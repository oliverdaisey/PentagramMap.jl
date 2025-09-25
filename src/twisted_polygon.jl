include("point.jl")
using Plots

const MonodromyMatrix = Matrix{BaseRing}

"""
    TwistedPolygon

A structure representing a twisted polygon in projective geometry. This consists of a sequence of vertices in homogeneous coordinates and a monodromy matrix that describes how the polygon "twists" as it wraps around.

Twisted polygons are more natural objects for the pentagram map than ordinary polygons, as they allow for a consistent definition of the map even when the polygon is not closed in the usual sense.
"""
struct TwistedPolygon
    vertices::Vector{Point}
    monodromy::MonodromyMatrix # 3x3 invertible matrix for projective transformation
end

# function TwistedPolygon(vertices::Vector{Point}, monodromy::MonodromyMatrix)
#     return TwistedPolygon(vertices, monodromy)
# end

"""
    twisted_polygon_from_vertices(vertices::Vector{Point}, monodromy::MonodromyMatrix) -> TwistedPolygon

Create a `TwistedPolygon` from a list of vertices and a monodromy matrix.
"""
function twisted_polygon_from_vertices(vertices::Vector{Point}, monodromy::MonodromyMatrix)
    if length(vertices) < 5
        error("A twisted polygon must have at least 5 vertices.")
    end
    if size(monodromy) != (3,3)
        error("Monodromy matrix must be 3x3.")
    end
    return TwistedPolygon(vertices, monodromy)
end

function Base.show(io::IO, tp::TwistedPolygon)
    print(io, "Twisted $(gonality(tp))-gon with monodromy ", monodromy(tp))
end

"""
    monodromy(tp::TwistedPolygon) -> MonodromyMatrix

Get the monodromy matrix of the twisted polygon.
"""
function monodromy(tp::TwistedPolygon)
    return tp.monodromy
end

"""
    vertices_modulo_monodromy(tp::TwistedPolygon) -> Vector{Point}

Get the vertices of the twisted polygon, considered modulo the action of the monodromy.
"""
function vertices_modulo_monodromy(tp::TwistedPolygon)
    return tp.vertices
end

"""
    gonality(tp::TwistedPolygon) -> Int

Get the gonality (number of vertices) of the twisted polygon.
"""
function gonality(tp::TwistedPolygon)
    return length(tp.vertices)
end

"""
    pentagram_map(tp::TwistedPolygon) -> TwistedPolygon

Compute the image of a twisted polygon under the pentagram map.
"""
function pentagram_map(tp::TwistedPolygon)
    n = gonality(tp)
    verts = vertices_modulo_monodromy(tp)
    newverts = Vector{Point}(undef, n)

    for i in 1:n
        v_im1 = verts[mod1(i-1, n)]
        v_i   = verts[i]
        v_ip1 = verts[mod1(i+1, n)]
        v_ip2 = verts[mod1(i+2, n)]

        l1 = line_through(v_im1, v_ip1)
        l2 = line_through(v_i, v_ip2)

        newverts[i] = normalize_point(intersection(l1, l2))
    end

    # monodromy update is subtle: in many formulations the new polygon inherits a conjugated monodromy. For now, we can keep the same matrix.
    return TwistedPolygon(newverts, monodromy(tp))
end

"""
    plot_twisted_polygon(tp::TwistedPolygon; cycles=1)

Plot the vertices of a twisted polygon in the affine plane.
You can specify how many monodromy cycles to draw.
"""
function plot_twisted_polygon(tp::TwistedPolygon; cycles=1)
    n = gonality(tp)
    verts = vertices_modulo_monodromy(tp)

    # Apply monodromy to extend the polygon
    points = Point[]
    for c in 0:cycles-1
        Mpow = monodromy(tp)^c
        for v in verts
            push!(points, Mpow * v)
        end
    end

    xy = map(affine_coords, points)
    xs = first.(xy)
    ys = last.(xy)

    plot(xs, ys, seriestype=:shape, lw=2, linecolor=:blue, fillalpha=0.1,
         marker=:circle, markersize=4, legend=false, aspect_ratio=1)
end

"""
    plot_pentagram_map(tp::TwistedPolygon, steps::Int; cycles=1)

Plot the twisted polygon after applying the pentagram map a specified number of times.
"""
function plot_pentagram_map(tp::TwistedPolygon, steps::Int; cycles=1)
    P = tp
    for i in 1:steps
        P = pentagram_map(P)
    end
    plot_twisted_polygon(P, cycles=cycles)
end

"""
    cross_ratio_coordinates(tp::TwistedPolygon)

Compute cross-ratio coordinates (x_i, y_i) for the twisted polygon `tp`.

Definitions used:
  p_i = intersection( (v_i, v_{i+2}), (v_{i+1}, v_{i+3}) )
  x_i = [ v_i, v_{i+1}, p_i, p_{i+1} ]
  y_i = [ v_{i+1}, v_{i+2}, p_i, p_{i+1} ]

Returns two vectors `(xs, ys)` of length n (same type as your numeric entries).
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