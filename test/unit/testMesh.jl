# include("../src/mesh/mesh2D.jl")
# include("../../src/mesh/mesh2D.jl")
# include("../../src/mesh/master2D.jl")
# include("../../src/mesh/mesh3D.jl")
# include("../../src/mesh/master3D.jl")
# include("../../src/mesh/compJacob.jl")

# mesh2d = Mesh2D( "square", 2 )

tol = 1e-13

# check boundaries
N = 5
M = 3
(p,t,bel) = makecube( M, N, M )

indchecks = [3, 1, 3, 1, 2, 2]
checkval  = [0.0, 1.0, 1.0, 0.0, 0.0, 1.0]

for jj in 1:size(bel,1)
  ic   = indchecks[ bel[jj,4] ]
  cval = checkval[ bel[jj,4] ]
  for ii in 1:3
    #@printf("%d %2.1f %d\n", jj, p[ bel[ jj, ii ], ic ], bel[jj,4] )
    @test p[ bel[ jj, ii ], ic ] == cval
  end
end

# check integration through volume integrals
for ii in [2,5,7,10], jj in [2,5,7,10], kk in [2,5,7,10], pp in 1:3

  mesh = Mesh3D( "cube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )

  compJacob!( mesh, master )

  vol = 0

  for el in 1:size(mesh.nodes,3)
    vol += sum( mesh.jcw[:,el] )
  end

  @test abs(vol - 1) < tol

end
