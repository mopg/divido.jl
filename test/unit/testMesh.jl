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

  @printf("%i %i %i %i\n",pp,ii,jj,kk)

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
