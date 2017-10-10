tol = 1e-13

# check boundaries
N = 5
M = 3
(p,t,bel) = luteos.makecube( M, N, M )

indchecks = [3, 1, 3, 1, 2, 2]
checkval  = [0.0, 1.0, 1.0, 0.0, 0.0, 1.0]

for jj in 1:size(bel,1)
  ic   = indchecks[ bel[jj,4] ]
  cval = checkval[ bel[jj,4] ]
  for ii in 1:3
    @test p[ bel[ jj, ii ], ic ] == cval
  end
end

# check connectivity
for ii in [2,5], jj in [2,5], kk in [2,5], pp in 1:3
  mesh   = Mesh3D( "cube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )
  for ff in 1:size(mesh.f,3)

    tr = mesh.f[ff,4]
    tl = mesh.f[ff,5]

    if tl > 0

      itrf = find( mesh.t2f[tr,1:4] .== ff )
      itlf = find( mesh.t2f[tl,1:4] .== ff )

      nodr = master.perm[ :,itrf,abs.(mesh.t2f[tr,itrf+4]) ]
      nodl = master.perm[ :,itlf,abs.(mesh.t2f[tl,itlf+4]) ]

      for qq in 1:length(nodr), dd in 1:3
        @test abs.(mesh.nodes[ nodr[qq], dd, tr ] -
                   mesh.nodes[ nodl[qq], dd, tl ]) < tol
      end

    else
      # check that nodes are actually on boundary they should be on
      itrf =   find( mesh.t2f[tr,1:4] .== ff )
      nodr =   master.perm[ :,itrf,abs.(mesh.t2f[tr,itrf+4]) ]
      nb   = - mesh.f[ff,end]
      for qq in 1:length(nodr)
        ic   = indchecks[ nb ]
        cval = checkval[  nb ]
        @test mesh.nodes[ nodr[qq], ic, tr ] == cval
      end

    end

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

# check integration through surface integrals
for ii in [5,10], jj in [2,10], kk in [2,10], pp in 1:3

  mesh = Mesh3D( "cube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )
  compJacob!( mesh, master )

  vol = 0

  for el in 1:size(mesh.nodes,3)

    # using divergence theorem: ∭ᵥ (∇⋅F) dV = ∯ₛ(F⋅n)dS
    # In this case F = 1/3 [x,y,z]

    voltemp = 0

    for qq in 1:4

      indF = mesh.t2f[el,qq]

      (ϕ2d,p2d,nods,normal,jcw2D) = luteos.compJacobFace( mesh, master, el, qq )

      volttemp = 0.0
      for ll in 1:size(p2d,1)
        volttemp += 1/3 * jcw2D[ll] * (normal[ll,:]' * p2d[ll,:])
      end
      voltemp += volttemp
    end

    vol += voltemp

  end

  @test abs(vol - 1) < tol

end

### Check skewed mesh
# check integration through volume integrals
for ii in [2,5,7,10], jj in [2,5,7,10], kk in [2,5,7,10], pp in 1:3

  mesh = Mesh3D( "skewcube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )

  compJacob!( mesh, master )

  vol = 0

  for el in 1:size(mesh.nodes,3)
    vol += sum( mesh.jcw[:,el] )
  end

  @test abs(vol - 1) < tol

end

# check integration through surface integrals
for ii in [5,10], jj in [2,10], kk in [2,10], pp in 1:3

  mesh = Mesh3D( "skewcube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )
  compJacob!( mesh, master )

  vol = 0

  for el in 1:size(mesh.nodes,3)

    # using divergence theorem: ∭ᵥ (∇⋅F) dV = ∯ₛ(F⋅n)dS
    # In this case F = 1/3 [x,y,z]

    voltemp = 0

    for qq in 1:4

      indF = mesh.t2f[el,qq]

      (ϕ2d,p2d,nods,normal,jcw2D) = luteos.compJacobFace( mesh, master, el, qq )

      volttemp = 0.0
      for ll in 1:size(p2d,1)
        volttemp += 1/3 * jcw2D[ll] * (normal[ll,:]' * p2d[ll,:])
      end
      voltemp += volttemp
    end

    vol += voltemp

  end

  @test abs(vol - 1) < tol

end

# check that basis functions and gauss points are correct
# for ii in [2,5,7,10], jj in [2,5,7,10], kk in [2,5,7,10], pp in 1:3
#
#   mesh = Mesh2D( "square", pp, N = ii )
#   master = Master2D( pp )
#
#   compJacob!( mesh, master )
#
#   barycent = fill( 0.0, size(master.gpts,1), 3 )
#   barycent[:,1]   = sum( master.gpts, 2 )
#   barycent[:,2:3] = master.gpts
#
#   for el in 1:size(mesh.nodes,3), ig in 1:size(master.gpts,1), dd in 1:2
#     cgpt = barycent[ig,1] * mesh.p[mesh.t[el,1],dd] +
#            barycent[ig,2] * mesh.p[mesh.t[el,2],dd] +
#            barycent[ig,3] * mesh.p[mesh.t[el,3],dd]
#     cgpt_ϕ = master.ϕ[:,ig]' * mesh.nodes[:,dd,el]
#     @test abs(cgpt_ϕ - cgpt) < tol
#   end
#
# end

# check integration through volume integrals using basis functions
# ∭ᵥ x dV
for ii in [2,5,7,10], jj in [2,5,7,10], kk in [2,5,7,10], pp in 1:3

  mesh = Mesh3D( "skewcube", pp, N = ii, M = jj, Q = kk )
  master = Master3D( pp )

  compJacob!( mesh, master )

  vol1 = 0
  vol2 = 0
  vol3 = 0

  for el in 1:size(mesh.nodes,3)
    vol1 += sum( diagm(mesh.jcw[:,el]) * (master.ϕ' * mesh.nodes[:,1,el] ) )
    vol2 += sum( diagm(mesh.jcw[:,el]) * (master.ϕ' * mesh.nodes[:,2,el] ) )
    vol3 += sum( diagm(mesh.jcw[:,el]) * (master.ϕ' * mesh.nodes[:,3,el] ) )
  end

  @test abs(vol1 - 1)   < tol
  @test abs(vol2 - 1)   < tol
  @test abs(vol3 - 0.5) < tol

end

# check that derivatives of basis functions are correct
for ii in [2,5,7,10], pp in 1:3
  mesh = Mesh2D( "square", pp, N = ii )
  master = Master2D( pp )

  compJacob!( mesh, master )

  for el in 1:size(mesh.nodes,3)
    ∇ϕc   = luteos.getderbfel( master, mesh.∂ξ∂x[:,el,:] )
    func  = mesh.nodes[:,1,el]
    res1  = ∇ϕc[:,:,1]' * func
    res2  = master.ϕ' * mesh.nodes[:,1,el].^0
    for cc in 1:length(res2)
      @test abs(res1[cc] - res2[cc]) < tol
    end

    func  = mesh.nodes[:,2,el]#0.5 * mesh.nodes[:,2,el].^2
    res1  = ∇ϕc[:,:,2]' * func
    res2  = master.ϕ' * mesh.nodes[:,2,el].^0
    for cc in 1:length(res2)
      @test abs(res1[cc] - res2[cc]) < tol
    end
  end

end
