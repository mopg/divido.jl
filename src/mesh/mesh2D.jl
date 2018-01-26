# ---------------------------------------------------------------------------- #
#
#   mesh2D.jl
#
#   Type for 2D meshes
#   Inherits from the abstract "mesh" type
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Mesh2D

Mesh2D type:
Type for 2D meshes consisting of triangles. It holds the nodes and connectivity
information.
"""
struct Mesh2D <: Mesh

  dim::Int64              # Dimension of the problem

  porder::Int64           # Polynomial order

  n::Int64                # Number of nodes

  p::Matrix{Float64}      # Nodal locations
  ploc::Matrix{Float64}   # Local nodal locations
  tloc::Matrix{Int64}     # Local triangles
  t::Matrix{Int64}        # Triangle - node connectivity
  t2f::Matrix{Int64}      # Triangle - face connectivity
  f::Matrix{Int64}        # Face - node/triangle connectivity
  fb::Matrix{Int64}       # Boundary face info
  nodes::Array{Float64,3} # Nodes on which solution is evaluated

end

"""
    Mesh2D( name::String, porder_::Int64; N = 5::Int64 )

Constructor for one of the default meshes. Currently "square" is implemented.
"""
function Mesh2D( name::String, porder::Porder; N = 5, M = N )

  porder_ = porder.p

  if name == "square"
    (p_, t_, bel_) = makesquare( N, M )
  elseif name == "single"
     p_   = [ 0.0 0.0;
              1.0 0.0;
              0.0 1.0 ]
     t_   = [ 1 2 3 ]
     bel_ = [ 2 3 1;
              3 1 2;
              1 2 3 ]
  elseif name[end-3:end] == ".su2"
    (p_, t_, bel_, tags_) = readSU2_2D( name )
  elseif name[end-3:end] == ".msh" # BAMG
    (p_, t_, bel_ ) = readBAMG( name )
  elseif name[end-4:end] == ".mesh" # FEFLOA
    (p_, t_, bel_ ) = readFEFLOA_2D( name )
  else
    error("Mesh2D: Unknown mesh type")
  end

  (f_, t2f_, nodes_, ploc_, tloc_, fb_) = genmesh( porder_, p_, t_, bel_ )

  n_ = size( p_, 1 )

  Mesh2D( 2, porder_, n_, p_, ploc_, tloc_, t_, t2f_, f_, fb_, nodes_)

end

"""
    Mesh2D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64}, bel::Matrix{Int64} )

Constructor for nodes locations (`p`), triangle connectivity (`t`), and
boundary element information (`bel`) given.
"""
function Mesh2D( porder::Porder, p_::Matrix{Float64}, t_::Matrix{Int64}, bel_::Matrix{Int64} )

  porder_ = porder.p

  (f_, t2f_, nodes_, ploc_, tloc_, fb_) = genmesh( porder_, p_, t_, bel_ )

  n_ = size( p_, 1 )

  Mesh2D( 2, porder_, n_, p_, ploc_, tloc_, t_, t2f_, f_, fb_, nodes_)

end

"""
    genmesh( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64}, bel::Matrix{Int64} )

Generates all necessary connectivity information given nodes locations (`p`),
triangle connectivity (`t`), and boundary element information (`bel`).
"""
function genmesh( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64}, bel_::Matrix{Int64} )

  (f_, t2f_, fb_)           = genFaces2D( t, bel_ )
  nodes_, ploc_, tloc_ = genNodes2D( porder, p, t )

  return f_, t2f_, nodes_, ploc_, tloc_, fb_

end

"""
    genFaces2D( t::Matrix{Int64}, bel::Matrix{Int64} )

Generates face connectivity given triangle connectivity (`t`), and
boundary element information (`bel`).
"""
function genFaces2D( t::Matrix{Int64}, bel::Matrix{Int64} )

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
  fb = fill( 0, size(bel,1), 2 )
  nfb = 1
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
    fb[nfb,1] = indf
    fb[nfb,2] = bel[indel[1],3]
    nfb += 1
  end

  # Look in fourth array
  indb = find( f[:,4] .== 0 )
  for jj = 1:length(indb)
    indf = indb[jj]

    # Find corresponding element in boundary element vector
    indel = find( all( bel[:,1:2] .== f[indf,1:2]', 2) )

    # Mark correct boundary
    f[indf,4] = - bel[indel[1],3]
    fb[nfb,1] = indf
    fb[nfb,2] = bel[indel[1],3]
    nfb += 1
  end

  return f, t2f, fb

end

"""
    genNodes2D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64} )

Generates the higher order nodes inside each element.
"""
function genNodes2D( porder::Int64, p::Matrix{Float64}, t::Matrix{Int64} )

  nt     = size(t,1)
  (plocal,tlocal) = genlocal( porder )

  npl    = size(plocal,1)
  nodes  = fill( 0.0, npl, 2, nt )

  # We use barycentric coordinates for this
  # for interpolation with barycentric coordinates at the point (xi,yi) with
  # barycentric coordinates (li1, li2, li3), the following holds:
  #   f(xi,yi) = li1 * f(x1,y1) + li2 * f(x2,y2) + li3 * f(x3,y3)

  for ii in 1:nt, jj in 1:npl
    for kk in 1:3
      nodes[jj,1,ii] += plocal[jj,kk] * p[ t[ii,kk], 1 ]
      nodes[jj,2,ii] += plocal[jj,kk] * p[ t[ii,kk], 2 ]
    end
  end

  return nodes, plocal, tlocal

end

"""
    genlocal( porder::Int64 )

Generates the barycentric coordinates for the higher order nodes inside
the master element (triangle).
"""
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

  # generate local triangles -- used for plotting
  tlocal = Matrix{Int64}( porder^2, 3 )
  if porder == 1
    tlocal[1,:]  = [1 2 3]
  elseif porder == 2
    tlocal[1,:]  = [1 2 4]
    tlocal[2,:]  = [2 5 4]
    tlocal[3,:]  = [2 3 5]
    tlocal[4,:]  = [4 5 6]
  elseif porder == 3
    tlocal[1,:]  = [1 8 7]
    tlocal[2,:]  = [8 10 7]
    tlocal[3,:]  = [8 9 10]
    tlocal[4,:]  = [9 4 10]
    tlocal[5,:]  = [9 2 4]
    tlocal[6,:]  = [7 10 6]
    tlocal[7,:]  = [10 5 6]
    tlocal[8,:]  = [10 4 5]
    tlocal[9,:]  = [6 5 3]
  elseif porder == 4
    tlocal[1,:]  = [1 10 9]
    tlocal[2,:]  = [10 13 9]
    tlocal[3,:]  = [10 11 13]
    tlocal[4,:]  = [11 14 13]
    tlocal[5,:]  = [11 12 14]
    tlocal[6,:]  = [12 4 14]
    tlocal[7,:]  = [12 2 4]
    tlocal[8,:]  = [9 13 8]
    tlocal[9,:]  = [13 15 8]
    tlocal[10,:] = [13 14 15]
    tlocal[11,:] = [14 5 15]
    tlocal[12,:] = [14 4 5]
    tlocal[13,:] = [8 15 7]
    tlocal[14,:] = [15 6 7]
    tlocal[15,:] = [15 5 6]
    tlocal[16,:] = [7 6 3]
  end

  return (plocal2, tlocal)

end

"""
    makesquare( n::Int64, m::Int64 )

Generates the nodes locations, triangle connectivity and boundary element
information for a square [0,1] x [0,1].
"""
function makesquare( n::Int64, m::Int64 )

  # boundary 1 is South
  # boundary 2 is East
  # boundary 3 is North
  # boundary 4 is West

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
  for ii in 1:n-1 # South
    bel[kk,:] = [ ii, ii + 1, 1 ]
    kk = kk + 1
  end
  for jj in 1:m-1 # East
    ii = n - 1
    bel[kk,:] = [ ii + 1 + (jj-1) * n, ii + 1 + jj * n, 2 ]
    kk = kk + 1
  end
  for ii in n-1:-1:1 # North
    jj = m - 1
    bel[kk,:] = [ ii + 1 + jj * n, ii + jj * n, 3 ]
    kk = kk + 1
  end
  for jj in m-1:-1:1 # East
    ii = 1
    bel[kk,:] = [ ii + jj * n, ii + (jj-1) * n, 4 ]
    kk = kk + 1
  end

  return p, t, bel

end
