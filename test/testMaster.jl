include("../src/mesh/master2D.jl")

tol = 1e-13

pmax = 4

for ii = 1:pmax

  @printf( " Test master element at %i\n", ii )

  master = Master2D( ii )

  println( "    1D")

  res1D = master.phi1D' * diagm(master.gwts1D) * master.phi1D
  sz    = size( res1D )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(res1D[jj,kk] - 1.0) < tol
    else
      @test abs(res1D[jj,kk]) < tol
    end
  end

  println( "    2D")

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
