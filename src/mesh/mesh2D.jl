using GroupSlices # TODO: KEEP TRACK OF DEPENDENCIES

type Mesh2D #<: Mesh    # Or should I just do union type??

  porder::Int64         # Polynomial order

  n::Int64              # Number of nodes

  p::Array{Float64}     # Nodal locations
  ploc::Array{Float64}  # Local nodal locations
  t::Array{Int64}       # Triangle - node connectivity
  t2f::Array{Int64}     # Triangle - face connectivity
  f::Array{Int64}       # Face - node/triangle connectivity
  nodes::Array{Float64} # Nodes on which solution is evaluated

  # Jacobian -- this doesn't seem like very efficient to store like this
  jcw::Array{Float64}
  ∂ξ∂x::Array{Float64}

end

# Constructor with name, i.e square
function Mesh2D( name::String, porder_::Int64; N = 5::Int64 )

  if name == "square"
    (p_, t_, bel_) = makesquare( N, N )
  else
    error("Unknown mesh type")
  end

  (f_, t2f_, nodes_, ploc_) = genmesh( porder_, p_, t_, bel_ )

  jcw_  = fill( 0.0, size(nodes_, 1), size(nodes_,3) )
  ∂ξ∂x_ = fill( 0.0, size(nodes_, 1), 4, size(nodes_,3) )

  n_ = size( p_, 1 )

  Mesh2D( porder_, n_, p_, ploc_, t_, t2f_, f_, nodes_, jcw_, ∂ξ∂x_ )

end

# Constructor with p, and t given

function Mesh2D( porder_::Int64, p_::Array{Float64}, t_::Array{Int64}, bel_::Array{Int64} )

  (f_, t2f_, nodes_, ploc_) = genmesh( porder_, p_, t_, bel_ )

  jcw_  = fill( 0.0, size(nodes_, 1), size(nodes_,3) )
  ∂ξ∂x_ = fill( 0.0, size(nodes_, 1), 4, size(nodes_,3) )

  n_ = size( p_, 1 )

  Mesh2D( porder_, n_, p_, ploc_, t_, t2f_, f_, nodes_, jcw_, ∂ξ∂x_ )

end

function genmesh( porder::Int64, p::Array{Float64}, t::Array{Int64}, bel_::Array{Int64} )

  (f_, t2f_)    = genFaces2D( t, bel_ )
  nodes_, ploc_ = genNodes2D( porder, p, t )

  return f_, t2f_, nodes_, ploc_

end

function genFaces2D( t::Array{Int64}, bel::Array{Int64} )

  nt = size(t,1)

  edges  = vcat( t[:,[2,3]], t[:,[3,1]], t[:,[1,2]] ) # This holds all edges (but multiple copies)
  triang = 1:nt
  triang = vcat( triang, triang, triang )

  ne = size(edges,1)

  tflip = Array{Int}( ne, 1 )

  for ii = 1:ne
    tflip[ii] = 1
    if edges[ii,1] > edges[ii,2]
      tflip[ii] = -1
      temp = edges[ii,1]
      edges[ii,1] = edges[ii,2]
      edges[ii,2] = temp
    end
  end

  # This index links to the first unique index in the array edges
  ix  = groupslices( edges, 1 )

  t2f = fill(0::Int64, nt, 3)

  # Find unique vector index
  # Need to find an index that links from edges to the unique edges (without gaps)
  jx = Array{Int64}( size(ix) )
  jx[1] = 0
  lind = 0 # maximum unique index
  mind = 0 # maximum index in edges that had unique index
  ne = 0  # number of unique elements

  indUni = fill( false, size(ix,1) ) # index of unique elements

  for ii = 1:length(ix)
    temp    = ix[ii] - lind - 1
    temp2   = ix[ii] - mind - 1
    jx[ii]  = temp
    mind    = max(mind,ix[ii])

    if temp2 >= 0 # Found a unique index
      lind = lind + 1
      ne  = ne + 1
      indUni[ii] = true
    end

    if temp2 < 0 # Need to link to the unique index
      jx[ii] = jx[ix[ii]]
    end

  end

  f = fill( 0::Int64, ne, 4 )

  f[:,1:2] = edges[indUni,:]

  # make t2f and complete f
  for jj = 1:length(tflip)

    ind  = triang[jj]
    indt = Int64(ceil(jj / nt))

    indf = ix[jj] - jx[jj]

    if tflip[jj] == 1
      f[indf,3] = ind
    else
      f[indf,4] = ind
    end

    t2f[ ind, indt ] = tflip[jj] * indf

  end

  ### Boundary Elements
  # include information from 'bel' into 'f': the index of the boundaries

  bel[:,1:2] = sort(bel[:,1:2],2) # sort such that we compare against f

  # Look in third array
  indb = find( f[:,3] .== 0 )
  for jj = 1:length(indb)
    indf = indb[jj]

    # Find corresponding element in boundary element vector
    indel = find( all( bel[:,1:2] .== f[indf,1:2]', 2) )

    # Ensure boundary element is oriented counter-clockwise
    #   flip orientation of face in t2f
    indt = f[indf,4]
    indfint = find( abs.(t2f[indt,:]) .== indf )
    t2f[ f[indf,4], indfint ] = - t2f[ f[indf,4], indfint ]
    #   flip nodes in f
    temp = f[indf,1]
    f[indf,1] = f[indf,2]
    f[indf,2] = temp
    f[indf,3] = f[indf,4]

    # Mark correct boundary
    f[indf,4] = - bel[indel[1],3]
  end

  # Look in fourth array
  indb = find( f[:,4] .== 0 )
  for jj = 1:length(indb)
    indf = indb[jj]

    # Find corresponding element in boundary element vector
    indel = find( all( bel[:,1:2] .== f[indf,1:2]', 2) )

    # Mark correct boundary
    f[indf,4] = - bel[indel[1],3]
  end

  return f, t2f

end

function genNodes2D( porder::Int64, p::Array{Float64}, t::Array{Int64} )

  nt     = size(t,1)
  # (plocal, tlocal) = genlocal( porder )
  plocal = genlocal( porder )

  npl    = size(plocal,1)
  nodes  = fill( 0.0, npl, 2, nt )

  # We use barycentric coordinates for this
  # for interpolation with barycentric coordinates at the point (xi,yi) with
  # barycentric coordinates (li1, li2, li3), the following holds:
  #   f(xi,yi) = li1 * f(x1,y1) + li2 * f(x2,y2) + li3 * f(x3,y3)

  for ii = 1:nt, jj = 1:npl
    for kk = 1:3
      nodes[jj,1,ii] += plocal[jj,kk] * p[ t[ii,kk], 1 ]
      nodes[jj,2,ii] += plocal[jj,kk] * p[ t[ii,kk], 2 ]
    end
  end

  return nodes, plocal

end

function genlocal( porder::Int64 )

  # Generate the barycentric coordinates
  # here these are implemented as (1 - x - y, x, y)

  n = Int64( round( 0.5 * (porder+1) * (porder+2) ) )

  plocal = Array{Float64}( n, 2 )
  kk = 1
  for jj = 1:porder+1
    for ii = 1:(porder + 2 - jj)
      plocal[kk,:] = [ (ii-1)/porder, (jj-1)/porder ]
      kk = kk + 1
    end
  end

  temp = 1 - sum(plocal,2)

  plocal = hcat( temp, plocal )

  # Switch to have corners first, then edges
  if porder > 1
    plocal2 = fill( 0.0, size(plocal) )

    #   corners
    plocal2[ 1:3, :] = plocal[ [1 porder+1 n], :]
    #   long edge
    kk = 2*porder + 1
    for jj in 1:(porder-1)
      plocal2[ (3+jj), :] = plocal[ kk, :]
      kk += porder - jj
    end
    #   west edge
    kk = n - 2
    for jj in 1:(porder-1)
      plocal2[ (2+porder+jj), :] = plocal[ kk, :]
      kk -= 2 + jj
    end
    #   south edge
    kk = 2
    for jj in 1:(porder-1)
      plocal2[ (1+2*porder+jj), :] = plocal[ kk, :]
      kk += 1
    end
    #   center
    kk = porder + 3
    pp = 3*(porder-1) + 4
    for jj in 1:(porder-2)
      for ii in 1:(porder-1-jj)
        plocal2[ pp, :] = plocal[ kk, :]
        pp += 1
        kk += 1
      end
      kk += 2
    end
  else
    plocal2 = plocal
  end

  return plocal2

end

function makesquare( n::Int64, m::Int64 )

  p   = Array{Float64}( n*m, 2 )
  t   = Array{Int64}( 2*(n-1)*(m-1), 3 )
  bel = Array{Int64}( 2*(n-1) + 2*(m-1), 3 )

  # nodes
  for ii = 1:n, jj = 1:m

    p[ii + (jj-1) * n, :] = [ 1/(n-1) * (ii - 1), 1/(m-1) * (jj - 1) ]

  end

  # triangles
  for ii = 1:(n-1), jj = 1:(m-1)

    # Lower triangles
    t[2*(ii + (jj-1) * (n-1)) - 1, :] = [ ii + (jj-1) * n,
              ii + 1 + (jj-1) * n,
              ii + 1 + jj     * n ]

    # Upper triangles
    t[2*(ii + (jj-1) * (n-1)),     :] = [ ii + (jj-1) * n,
              ii + 1 + jj     * n,
              ii     + jj     * n ]

  end

  # boundary elements
  kk = 1
  for ii = 1:n-1 # South
    bel[kk,:] = [ ii, ii + 1, 1 ]
    kk = kk + 1
  end
  for jj = 1:m-1 # East
    ii = n - 1
    bel[kk,:] = [ ii + 1 + (jj-1) * n, ii + 1 + jj * n, 2 ]
    kk = kk + 1
  end
  for ii = n-1:-1:1 # North
    jj = m - 1
    bel[kk,:] = [ ii + 1 + jj * n, ii + jj * n, 3 ]
    kk = kk + 1
  end
  for jj = m-1:-1:1 # East
    ii = 1
    bel[kk,:] = [ ii + jj * n, ii + (jj-1) * n, 4 ]
    kk = kk + 1
  end

  return p, t, bel

end
