# ---------------------------------------------------------------------------- #
#
#   mesh3D.jl
#
#   Type for 3D meshes
#   Inherits from the abstract "mesh" type
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Mesh3D

Mesh3D type:
Type for 3D meshes consisting of tetrahedrons. It holds the nodes and
connectivity information.
"""
type Mesh3D <: Mesh

  dim::Int64            # Dimension of the problem

  porder::Int64         # Polynomial order

  n::Int64              # Number of nodes

  p::Matrix{Float64}     # Nodal locations
  ploc::Matrix{Float64}  # Local nodal locations
  tloc::Matrix{Int64}    # Local tets
  t::Matrix{Int64}       # Tet - node connectivity
  t2f::Matrix{Int64}     # Tet - face connectivity
  f::Matrix{Int64}       # Face - node/tet connectivity
  nodes::Array{Float64,3} # Nodes on which solution is evaluated
  trorder::Matrix{Int64} # Different ordering indices for faces

  # Jacobian -- this doesn't seem like very efficient to store like this
  jcw::Matrix{Float64}
  ∂ξ∂x::Array{Float64,3}

end

"""
    Mesh3D( name::String, porder_::Int64; N = 5::Int64, M = N, Q = N )

Constructor for one of the default meshes. Currently "cube" is implemented.
"""
function Mesh3D( name::String, porder_::Int64; N = 5::Int64, M = N, Q = N )

  setup()

  if name == "cube"
    (p_, t_, bel_) = makecube( N, M, Q )
  else
    error("Unknown mesh type")
  end

  if porder_ > 3
    error("P>3 not implemented for 3D")
  end

  (f_, t2f_, nodes_, ploc_, tloc_, trorder_) = genmesh3D( porder_, p_, t_, bel_ )

  jcw_  = fill( 0.0, size(nodes_, 1), size(nodes_,3) )
  ∂ξ∂x_ = fill( 0.0, size(nodes_, 1), 9, size(nodes_,3) )

  n_ = size( p_, 1 )

  Mesh3D( 3, porder_, n_, p_, ploc_, tloc_, t_, t2f_, f_, nodes_, trorder_, jcw_, ∂ξ∂x_ )

end

"""
    genmesh3D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64}, bel::Matrix{Int64} )

Generates all necessary connectivity information given nodes locations (`p`),
tetrahedron connectivity (`t`), and boundary element information (`bel`).
"""
function genmesh3D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64}, bel::Matrix{Int64} )

  (f_, t2f_, trorder_) = genFaces3D( t, bel )
  nodes_, ploc_, tloc_ = genNodes3D( porder, p, t )

  return f_, t2f_, nodes_, ploc_, tloc_, trorder_

end

"""
    genFaces3D( t::Matrix{Int64}, bel::Matrix{Int64} )

Generates face connectivity given tetrahedron connectivity (`t`), and
boundary element information (`bel`).
"""
function genFaces3D( t::Matrix{Int64}, bel::Matrix{Int64} )

  nt = size(t,1)

  # 2 3 4 - 1 4 3  - 1 2 4 - 1 3 2 : First ind is the vertices of face which does not included current vertices
  faces = vcat( t[:,[2,3,4]], t[:,[1,4,3]], t[:,[1,2,4]], t[:,[1,3,2]] ) # This holds all faces (but multiple copies)
  tets  = 1:nt
  tets  = vcat( tets, tets, tets, tets )

  nf = size(faces,1)

  orderind = [1 2 3;
              1 3 2;
              2 1 3;
              2 3 1;
              3 1 2;
              3 2 1]
  orderflip = [1, -1, -1, 1, 1, -1]

  tflip     = fill( 1::Int64, nf )
  typeorder = fill( 1::Int64, nf )

  # sort in ascending order
  for ii in 1:nf
    indorder      = sortperm(  faces[ii,:] )
    ind           = find( all( indorder' .== orderind, 2 ) )
    tflip[ii]     = orderflip[ ind[1] ]
    typeorder[ii] = ind[1]
    faces[ii,:]   = faces[ii, orderind[ind[1],:] ]
  end

  # This index links to the first unique index in the array faces
  ix  = groupslices( faces, 1 )

  t2f     = fill(0::Int64, nt, 8)
  tforder = fill(0::Int64, nt, 3)

  # Find unique vector index
  # Need to find an index that links from faces to the unique faces (without gaps)
  jx = Array{Int64}( size(ix) )
  jx[1] = 0
  lind  = 0 # maximum unique index
  mind  = 0 # maximum index in faces that had unique index
  nf    = 0  # number of unique elements

  indUni = fill( false, size(ix,1) ) # index of unique elements

  for ii = 1:length(ix)
    temp    = ix[ii] - lind - 1
    temp2   = ix[ii] - mind - 1
    jx[ii]  = temp
    mind    = max(mind,ix[ii])

    if temp2 >= 0 # Found a unique index
      lind = lind + 1
      nf  += 1
      indUni[ii] = true
    end

    if temp2 < 0 # Need to link to the unique index
      jx[ii] = jx[ix[ii]]
    end

  end

  f = fill( 0::Int64, nf, 5 )

  f[:,1:3] = faces[indUni,:]

  # make t2f and complete f
  for jj = 1:length(tflip)

    ind  = tets[jj]
    indt = Int64(ceil(jj / nt))

    indf = ix[jj] - jx[jj]

    if tflip[jj] == 1
      f[indf,4] = ind
    else
      f[indf,5] = ind
    end

    #t2f[ ind, indt ] = tflip[jj] * indf # NOTE: JUST A MINUS SIGN IS NOT SUFFICIENT ANYMORE, NEED TO KNOW THE ORDERING EXACTLY (I.E. ONE OF SIX DIFFERENT ORDERINGS)
    t2f[ ind, indt ]     = indf
    t2f[ ind, indt + 4 ] = tflip[jj] * typeorder[jj] # Sign indicates whether outward normal or not (outward: pos)

  end

  ### Boundary Elements
  # include information from 'bel' into 'f': the index of the boundaries

  bel[:,1:3] = sort(bel[:,1:3],2) # sort such that we compare against f

  # Look in fourth column of array
  indb = find( f[:,4] .== 0 )
  for jj = 1:length(indb)
    indf = indb[jj]

    # Find corresponding element in boundary element vector
    indel = find( all( bel[:,1:3] .== f[indf,1:3]', 2) )

    # Ensure boundary element is oriented counter-clockwise
    #   flip orientation of face in t2f
    indt     = f[indf,5]
    indfint  = find( t2f[indt,1:4] .== indf )
    orderorg = orderind[ abs.(t2f[indt,4+indfint]), : ]
    # flip nodes to get outward normal
    temp        = orderorg[3]
    orderorg[3] = orderorg[2]
    orderorg[2] = temp
    indorder = find( all( orderorg .== orderind, 2 ) )
    t2f[indt,4+indfint] = indorder
    #   flip nodes in f
    temp = f[indf,2]
    f[indf,2] = f[indf,3]
    f[indf,3] = temp
    f[indf,4] = f[indf,5]

    # Mark correct boundary
    f[indf,5] = - bel[indel[1],4]
  end

  # Look in fourth array
  indb = find( f[:,5] .== 0 )
  for jj = 1:length(indb)
    indf = indb[jj]

    # Find corresponding element in boundary element vector
    indel = find( all( bel[:,1:3] .== f[indf,1:3]', 2) )

    # Mark correct boundary
    f[indf,5] = - bel[indel[1],4]
  end

  return f, t2f, orderind

end

"""
    genNodes3D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64} )

Generates the higher order nodes inside each element.
"""
function genNodes3D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64} )

  nt     = size(t,1)
  # (plocal, tlocal) = genlocal( porder )
  (plocal,tlocal) = genlocal3D( porder )

  npl    = size(plocal,1)
  nodes  = fill( 0.0, npl, 3, nt )

  # We use barycentric coordinates for this
  # for interpolation with barycentric coordinates at the point (xi,yi,zi) with
  # barycentric coordinates (li1, li2, li3, li3), the following holds:
  #   f(xi,yi) = li1 * f(x1,y1) + li2 * f(x2,y2) + li3 * f(x3,y3) + li4 * f(x4,y4)

  for ii in 1:nt, jj in 1:npl
    for kk in 1:4
      nodes[jj,1,ii] += plocal[jj,kk] * p[ t[ii,kk], 1 ]
      nodes[jj,2,ii] += plocal[jj,kk] * p[ t[ii,kk], 2 ]
      nodes[jj,3,ii] += plocal[jj,kk] * p[ t[ii,kk], 3 ]
    end
  end

  return nodes, plocal, tlocal

end

"""
    genlocal3D( porder::Int64 )

Generates the barycentric coordinates for the higher order nodes inside
the master element (tetrahedron).
"""
function genlocal3D( porder::Int64 )

  # The barycentric coordinates are implemented as (1 - x - y - z, x, y, z)

  n   = Int64( round( 1/6 * (porder+1) * (porder+2) * (porder+3) ) )
  n2d = Int64( round( 0.5 * (porder+1) * (porder+2) ) )

  plocal = Array{Float64}( n, 3 )
  qq = 1
  for kk = 1:porder+1
    for jj = 1:porder+2-kk
      for ii = 1:(porder + 3 - jj - kk)
        plocal[qq,:] = [ (ii-1)/porder, (jj-1)/porder, (kk-1)/porder ]
        qq += 1
      end
    end
  end

  temp = 1 - sum(plocal,2)

  plocal = hcat( temp, plocal )

  # Switch to have corners first, then edges
  plocal2 = fill( 0.0, size(plocal) )
  #COORD(3, s) = {0, 1, 0, 0, 0 , 0 , O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, T3, O3, 0 , O3, O3};
  #COORD(3, t) = {0, 0, 1, 0, T3, O3, 0 , 0 , O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, O3, 0 , O3};
  #COORD(3, u) = {0, 0, 0, 1, O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, T3, 0 , 0 , O3, O3, O3, 0 };
  indmap = [1,2,3,4]
  if porder == 2
    indmap = [1,porder+1,n2d,n,9,8,5,4,7,2]
  elseif porder == 3
    indmap = [1,porder+1,n2d,n,16,19,18,13,7,9,8,5,11,17,2,3,15,14,12,6]
  end
  plocal2 = plocal[ indmap, :]

  # generate local tetrahedra -- used for plotting
  if porder == 1
    tlocal      = fill( 0, 1, 4)
    tlocal[1,:] = [1 2 3 4]
  elseif porder == 2
    tlocal = [10     6     2     7
               4     9     5     6
               8     5     9    10
               3     5     8     7
               6     9     5    10
               8     5    10     7
               6    10     5     7
               8     9     1    10]
  elseif porder == 3
    tlocal = [2    16     8     9
             18     5    17    20
             20     5    17    10
             18    15    12    20
             18    19    15    20
             20    17     9    10
              3     5    11    10
             17     7    19     8
             18    17    19    20
             11     5    20    10
             16    17     8     9
              6     4    14     7
             11     5    18    20
             11    18    12    20
             20    17    16     9
             18     5     6    17
             15    20    19    16
             16    17    19     8
             20    17    19    16
             18    14    13    19
             19     6    14     7
             19    17     6     7
             18    17     6    19
             18     6    14    19
             18    13    15    19
             18    13    12    15
              1    12    13    15]
  end

  return (plocal2, tlocal)

end

"""
    makecube( n::Int64, m::Int64, q::Int64 )

Generates the nodes locations, tetrahedron connectivity and boundary element
information for a cube [0,1] x [0,1] x [0,1].
"""
function makecube( n::Int64, m::Int64, q::Int64 )

  # boundary 1 is z=0
  # boundary 2 is y=1
  # boundary 3 is z=1
  # boundary 4 is y=0
  # boundary 5 is x=0
  # boundary 6 is x=1

  p   = Array{Float64}( n*m*q, 3 )
  hex = Array{Int64}( (n-1)*(m-1)*(q-1), 8 ) # Map to ordering per hex cell
  t   = Array{Int64}( 6*(n-1)*(m-1)*(q-1), 4 )
  bel = Array{Int64}( 4*(n-1)*(q-1) + 4*(m-1)*(q-1) + 4*(m-1)*(n-1), 4 )

  # nodes
  for ii = 1:n, jj = 1:m, kk = 1:q
    p[ii + (jj-1) * n + (kk-1) * n * m, :] = [ 1/(n-1) * (ii - 1),
      1/(m-1) * (jj - 1), 1/(q-1) * (kk - 1) ]
  end
  # hex nodes
  kk     = 1
  hexind = 1
  for qq = 1:(q-1)
    for jj = 1:(m-1)
      for ii = 1:(n-1)

        hex[ kk, : ] = [hexind,     hexind+1,     hexind+n,     hexind+n+1,
                        hexind+n*m, hexind+n*m+1, hexind+n+n*m, hexind+n+1+n*m]

        hexind += 1
        kk += 1

      end
      hexind += 1
    end
    hexind += n
  end

  #         y
  #  3----------4
  #  |\     ^   |\
  #  | \    |   | \
  #  |  \   |   |  \
  #  |   7------+---8
  #  |   |  +-- |-- | -> x
  #  1---+---\--2   |
  #   \  |    \  \  |
  #    \ |     \  \ |
  #     \|      z  \|
  #      5----------6


  # tet index within hex
  tethexind = [1 2 3 5;
               2 5 4 3;
               7 3 4 5;
               2 5 6 4;
               5 7 6 4;
               8 7 4 6]

  # tetrahedra
  ntet_level = (n-1) * (m-1)
  for ii = 1:(n-1), jj = 1:(m-1), kk = 1:(q-1)

    hexind   = (ii-1) + (jj-1) * (n-1) + (kk-1) * (n-1) * (m-1)

    for pp in 1:6
      t[ 6*hexind + pp, :] = hex[ hexind+1, tethexind[pp,:] ]
    end

  end

  # triangle index for boundary
  indyzpos  = [2 4 6;
               4 8 6]
  indyzneg  = [1 5 3;
               5 7 3]
  indxzneg  = [6 5 2;
               5 1 2]
  indxzpos  = [8 4 7;
               4 3 7]
  indxypos  = [5 6 7;
               6 8 7]
  indxyneg  = [4 2 3;
               2 1 3]

  # boundary elements
  kk = 1
  for ii in 1:n-1, jj in 1:m-1 # xy neg
    indhex = (ii-1) + (jj-1)*(n-1)
    bel[kk,:] = [ hex[ indhex+1, indxyneg[1,:] ]'  1 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indxyneg[2,:] ]'  1 ]
    kk = kk + 1
  end
  for jj in 1:m-1, qq in 1:q-1 # yz pos
    indhex = (n-2) + (jj-1)*(n-1) + (qq-1)*(n-1)*(m-1)
    bel[kk,:] = [ hex[ indhex+1, indyzpos[1,:] ]'  2 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indyzpos[2,:] ]'  2 ]
    kk = kk + 1
  end
  for ii in 1:n-1, jj in 1:m-1 # xy pos
    indhex = (ii-1) + (jj-1)*(n-1) + (q-2)*(n-1)*(m-1)
    bel[kk,:] = [ hex[ indhex+1, indxypos[1,:] ]'  3 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indxypos[2,:] ]'  3 ]
    kk = kk + 1
  end
  for jj in 1:m-1, qq in 1:q-1 # yz neg
    indhex = (jj-1)*(n-1) + (qq-1)*(n-1)*(m-1)
    bel[kk,:] = [ hex[ indhex+1, indyzneg[1,:] ]'  4 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indyzneg[2,:] ]'  4 ]
    kk = kk + 1
  end
  for ii in 1:n-1, qq in 1:q-1 # xz neg
    indhex = (ii-1) + (qq-1)*(n-1)*(m-1)
    bel[kk,:] = [ hex[ indhex+1, indxzneg[1,:] ]'  5 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indxzneg[2,:] ]'  5 ]
    kk = kk + 1
  end
  for ii in 1:n-1, qq in 1:q-1 # xz pos
    indhex = (ii-1) + (m-2)*(n-1) + (qq-1)*(n-1)*(m-1)
    bel[kk,:] = [ hex[ indhex+1, indxzpos[1,:] ]'  6 ]
    kk = kk + 1
    bel[kk,:] = [ hex[ indhex+1, indxzpos[2,:] ]'  6 ]
    kk = kk + 1
  end

  return p, t, bel

end
