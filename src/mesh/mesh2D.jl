using GroupSlices # TODO: KEEP TRACK OF DEPENDENCIES

type Mesh2D #<: Mesh # Or should I just do union type??

  porder::Int64       # Polynomial order

  n::Int64            # Number of nodes

  p::Array{Float64}   # Nodal locations
  t::Array{Int64}     # Triangle - node connectivity
  t2f::Array{Int64}   # Triangle - face connectivity
  f::Array{Int64}     # Face - node/triangle connectivity
  nodes::Array{Int64} # Nodes on which solution is evaluated

end

# Constructor with name, i.e square
function Mesh2D( name::String, porder_::Int64 )

  if name == "square"
    (p_, t_) = makesquare( 5 )
  else
    error("Unknown mesh type")
  end

  (f_, t2f_, nodes_) = genmesh( porder_, p_, t_ )

  n_ = size( p_, 1 )

  Mesh2D( porder_, n_, p_, t_, t2f_, f_, nodes_ )

end

# Constructor with p, and t given

function Mesh2D( porder_::Int64, p_::Array{Float64}, t_::Array{Int64} )

  (f_, t2f_, nodes_) = genmesh( porder_, p_, t_ )

  n_ = size( p_, 1 )

  Mesh2D( porder_, n_, p_, t_, t2f_, f_, nodes_ )

end

function genmesh( porder::Int64, p::Array{Float64}, t::Array{Int64} )

  f_     = [0]
  t2f_   = [0]
  nodes_ = [0]
  (f_, t2f_) = genFaces2D( t ) # TODO: Boundaries! (have that be an input)
  # nodes_     = genNodes2D( porder, p, t )
  nodes_ = [0]
  return f_, t2f_, nodes_

end

function genFaces2D( t::Array{Int64} )

  nt = size(t,1)

  edges  = vcat( t[:,[2,3]], t[:,[3,1]], t[:,[1,2]] )
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

  ix     = groupslices( edges, 1 )

  t2f = fill(0::Int64, nt, 3)

  # Find unique vector index
  jx = Array{Int64}( size(ix) )
  jx[1] = 0
  lind = 0
  mind = 0

  ne = 0

  indUni = fill( false, size(ix,1) )

  for ii = 1:length(ix)
    temp = ix[ii] - lind - 1;temp2 = ix[ii] - mind - 1
    jx[ii] = temp
    mind = max(mind,ix[ii])
    if temp >= 0
      lind = lind + 1
      ne  = ne + 1
      indUni[ii] = true
    end
    if temp2 < 0
      jx[ii] = jx[ix[ii]]
    end
  end

  println(ix)
  println(tflip)
  println(ix - jx)

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

  return f, t2f

end

function makesquare( a )

  n = 3
  m = 3

  p = Array{Float64}( n*m, 2 )
  t = Array{Int64}( 2*(n-1)*(m-1), 3 )

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

  return p, t

end
