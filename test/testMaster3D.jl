# include("../src/mesh/master3D.jl")

tol = 1e-13

pmax = 3

for ii = 1:pmax

  @printf( " Test 3D master element at %i\n", ii )

  master = Master3D( ii )

  println( "    3D")

  res = master.phi' * diagm(master.gwts) * master.phi
  sz  = size( res )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(res[jj,kk] - 1.0) < tol
    else
      @test abs(res[jj,kk]) < tol
    end
  end

end
