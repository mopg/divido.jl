type Mesh3D #<: Mesh

  dim::Int64            # Dimension of the problem

  porder::Int64         # Polynomial order

  n::Int64              # Number of nodes

  p::Array{Float64}     # Nodal locations
  ploc::Array{Float64}  # Local nodal locations
  tloc::Array{Int64}    # Local tets
  t::Array{Int64}       # Tet - node connectivity
  t2f::Array{Int64}     # Tet - face connectivity
  f::Array{Int64}       # Face - node/tet connectivity
  nodes::Array{Float64} # Nodes on which solution is evaluated

  # Jacobian -- this doesn't seem like very efficient to store like this
  jcw::Array{Float64}
  ∂ξ∂x::Array{Float64}

end

# Constructor with name, i.e cube
function Mesh3D( name::String, porder_::Int64; N = 5::Int64, M = N, Q = N )

  if name == "cube"
    (p_, t_, bel_) = makecube( N, M, Q )
  else
    error("Unknown mesh type")
  end

  println(bel_)

  f_     = fill(0.0,2,2)
  t2f_   = fill(0.0,2,2)
  nodes_ = fill(0.0,2,2)
  ploc_  = fill(0.0,2,2)
  tloc_  = fill(0.0,2,2)
  # (f_, t2f_, nodes_, ploc_, tloc_) = genmesh3D( porder_, p_, t_, bel_ )

  jcw_  = fill( 0.0, size(nodes_, 1), size(nodes_,3) )
  ∂ξ∂x_ = fill( 0.0, size(nodes_, 1), 9, size(nodes_,3) )

  n_ = size( p_, 1 )

  Mesh3D( 3, porder_, n_, p_, ploc_, tloc_, t_, t2f_, f_, nodes_, jcw_, ∂ξ∂x_ )

end

function genmesh3D( porder::Int64, p::Array{Float64}, t::Array{Int64}, bel::Array{Int64} )

  (f_, t2f_)           = genFaces3D( t, bel )
  nodes_, ploc_, tloc_ = genNodes3D( porder, p, t )

  return f_, t2f_, nodes_, ploc_, tloc_

end

function makecube( n::Int64, m::Int64, q::Int64 )

  # boundary 1 is South
  # boundary 2 is East
  # boundary 3 is North
  # boundary 4 is West
  # boundary 5 is South
  # boundary 6 is North

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

  println(hex)

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
