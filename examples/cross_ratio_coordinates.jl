verts = [
    Rational{Int}[1,0,1],
    Rational{Int}[1,1,1],
    Rational{Int}[0,1,1],
    Rational{Int}[-1,0,1],
    Rational{Int}[0,-1,1]
]
M = Matrix{Rational{Int}}(I, 3, 3)
tp = TwistedPolygon(verts, M)
xs, ys = cross_ratio_coordinates(tp)
println("x = ", xs)
println("y = ", ys)

tp2 = pentagram_map(tp)
xs2, ys2 = cross_ratio_coordinates(tp2)
println("After one pentagram step:")
println("x = ", xs2)
println("y = ", ys2)