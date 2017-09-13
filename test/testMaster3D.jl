# include("../src/mesh/master3D.jl")

tol = 1e-13

pmax = 3

println( "Legendre basis functions")

for ii = 1:pmax

  @printf( "   Test 3D master element at %i\n", ii )

  master = Master3D( ii; typeb = "leg" )

  println( "      3D")

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

# TODO: Lagrangian
