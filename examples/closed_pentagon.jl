include("../src/twisted_polygon.jl")
 
 # Define a twisted pentagon with a specific monodromy

verts = [
           BaseRing[1, 0, 1],
           BaseRing[2, 2, 1],
           BaseRing[-2, 1, 1],
           BaseRing[-1, -1, 1],
           BaseRing[0, -1, 1]
       ]
M = Matrix{BaseRing}([
           0 -1 0;
           1 0 0;
           0 0 1
       ])
M = Matrix{BaseRing}([
           1 0 0;
           0 1 0;
           0 0 1
       ])
tp = twisted_polygon_from_vertices(verts, M)
plot_twisted_polygon(tp)