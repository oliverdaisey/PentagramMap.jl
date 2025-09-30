using PentagramMap

verts = [
    Rational{Int}[1,0,1],
    Rational{Int}[1,1,1],
    Rational{Int}[0,1,1],
    Rational{Int}[-1,0,1],
    Rational{Int}[-2,-2,1],
    Rational{Int}[0,-1,1]
]
M = [1 0 0;
     0 1 0;
     0 0 1]
tp = PentagramMap.TwistedPolygon(verts, M)
plot_twisted_polygon(tp)
xs, ys = cross_ratio_coordinates(tp)
println("x = ", xs)
println("y = ", ys)

tp2 = pentagram_map(tp)
xs2, ys2 = cross_ratio_coordinates(tp2)
println("After one pentagram step:")
println("x = ", xs2)
println("y = ", ys2)

tp3 = pentagram_map(tp2)
xs3, ys3 = cross_ratio_coordinates(tp3)
println("After two pentagram steps:")
println("x = ", xs3)
println("y = ", ys3)

println("Using cross ratio coordinates formula:")
cxs2, cys2 = pentagram_map(xs, ys)

println("Difference in x: ", xs2 - cxs2)
println("Difference in y: ", ys2 - cys2)

# we check that the casimirs are preserved

O_n, E_n = casimirs(xs, ys)
O_n2, E_n2 = casimirs(xs2, ys2)

println("Casimirs before pentagram step: O_n = $O_n, E_n = $E_n")
println("Casimirs after pentagram step: O_n = $O_n2, E_n = $E_n2")