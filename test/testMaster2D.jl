# include("../src/mesh/master2D.jl")
# include("../src/mesh/mesh2D.jl")

tol = 1e-13

pmax = 4

### Legendre basis functions
println( "Lagrangian basis functions")
for ii = 1:pmax

  ploc   = genlocal(ii)
  ploc1D = fill( 0.0, ii+1, 1 )
  for jj = 1:ii+1
    ploc1D[jj] = (jj-1)/ii
  end
  temp = ploc1D[2:ii]
  ploc1D[1:2] = ploc1D[ [1,ii+1] ]
  ploc1D[3:ii+1] = temp

  @printf( "  Test 2D master element at %i\n", ii )

  master = Master2D( ii; typeb = "lag" )

  println( "      1D")

  phi1D, dphi1D = basisFuncLineLag( Val{ii}, ploc1D )
  sz    = size( phi1D )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(phi1D[jj,kk] - 1.0) < tol
    else
      @test abs(phi1D[jj,kk]) < tol
    end
  end

  println( "      2D")

  phi, dphi = basisFuncTriangleLag( Val{ii}, ploc[:,2], ploc[:,3] )

  sz  = size( phi )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(phi[jj,kk] - 1.0) < tol
    else
      @test abs(phi[jj,kk]) < tol
    end
  end

end

### Legendre basis functions
println( "Legendre basis functions")
for ii = 1:pmax

  @printf( "  Test 2D master element at %i\n", ii )

  master = Master2D( ii; typeb = "leg" )

  println( "      1D")

  res1D = master.ϕ1D' * diagm(master.gwts1D) * master.ϕ1D
  sz    = size( res1D )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(res1D[jj,kk] - 1.0) < tol
    else
      @test abs(res1D[jj,kk]) < tol
    end
  end

  println( "      2D")

  res = master.ϕ' * diagm(master.gwts) * master.ϕ
  sz  = size( res )
  for jj = 1:sz[1], kk = 1:sz[2]
    if jj == kk
      @test abs(res[jj,kk] - 1.0) < tol
    else
      @test abs(res[jj,kk]) < tol
    end
  end

end
