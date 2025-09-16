using StaticArrays
using Plots

# cross product in homogeneous coords
cross3(a::SVector{3,T}, b::SVector{3,T}) where T = SVector(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
)

r"""
    pentagram_step(polygon::Vector{SVector{3,Float64}}) -> Vector{SVector{3,Float64}}

Perform one step of the pentagram map on a polygon given in homogeneous coordinates.
"""
function pentagram_step(polygon::Vector{SVector{3,T}}) where T
    n = length(polygon)
    newpoly = Vector{SVector{3,T}}(undef, n)
    for i in 1:n
        p1 = polygon[i]
        p2 = polygon[mod1(i+2,n)]
        p3 = polygon[mod1(i+1,n)]
        p4 = polygon[mod1(i+3,n)]
        l1 = cross3(p1,p2)
        l2 = cross3(p3,p4)
        newpt = cross3(l1,l2)
        newpoly[i] = newpt
    end
    return newpoly
end

function plot_polygon(polygon::Vector{SVector{3,T}}) where T
    n = length(polygon)
    xs = Float64[]
    ys = Float64[]
    for i in 1:n
        p = polygon[i]
        push!(xs, p[1]/p[3])
        push!(ys, p[2]/p[3])
    end
    push!(xs, xs[1])
    push!(ys, ys[1])
    plot(xs, ys, seriestype=:shape, aspect_ratio=:equal, legend=false)
end

function plot_pentagram_map(polygon::Vector{SVector{3,T}}, steps::Int) where T

    P = polygon
    for i in 0:steps
        P = pentagram_step(P)
        println(P)
    end
    plot_polygon(P)
end