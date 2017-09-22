# include("../src/mesh/mesh2D.jl")
include("../../src/mesh/mesh3D.jl")

# mesh2d = Mesh2D( "square", 2 )

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
