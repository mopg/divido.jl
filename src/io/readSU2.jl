# ---------------------------------------------------------------------------- #
#
#   readSU2.jl
#
#   Read SU2 meshes into luteos
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    readSU2_2D( flname::String )

Reads in 2D SU2 mesh.
"""
function readSU2_2D( flname::String )

  dim = 2

  # Open file
  fid = open( flname, "r" )

  # Check dimension
  x = "%"
  while x[1] == '%'
      x = readline(fid)
  end
  xx   = split(x)
  ndim = parse( Int64, xx[end] )
  if ndim != dim
      error(" Input mesh error: Input mesh is wrong dimension, should be 2D")
  end

  # Read number of elements and triangle connectivity
  x = "%"
  while x[1] == '%'
      x = readline(fid)
  end
  xx    = split( x )
  nelem = parse( Int64, xx[end] )
  t = fill( 0, nelem, 3 )
  for il in 1:nelem #x[1] != '%' || x[1] != 'N'
      x   = readline(fid)
      xx  = split(x)
      ind = parse(Int64,xx[end]) + 1

      t[ind,1] = parse(Int64,xx[2]) + 1
      t[ind,2] = parse(Int64,xx[3]) + 1
      t[ind,3] = parse(Int64,xx[4]) + 1
  end
  if parse(Int64,xx[1]) != 5
      error(" Input mesh error: Only triangles are supported for 2D")
  end

  # Read node locations
  x = "%"
  while x[1] == '%'
      x = readline(fid)
  end
  xx = split( x )
  np = parse( Int64, xx[end] )
  p  = fill( 0.0, np, 2 )
  for il in 1:np
      x   = readline(fid)
      xx  = split(x)
      ind = parse(Int64,xx[end]) + 1

      p[ind,1] = parse(Float64,xx[1])
      p[ind,2] = parse(Float64,xx[2])
  end

  # Read boundary elements
  x = "%"
  while x[1] == '%'
      x = readline(fid)
  end
  xx = split( x )
  nm = parse( Int64, xx[end] )
  nnm = fill( 0, nm )
  sm = fill( "tag", nm )
  belm = fill( fill(0, 0, 0), nm )
  nelm = fill( 0, nm )

  for jjb in 1:nm
      x         = readline( fid )
      xx        = split( x )
      sm[jjb]   = xx[end]

      x         = readline( fid )
      xx        = split( x )
      nel       = parse( Int64, xx[end] )
      nelm[jjb] = nel
      belm[jjb] = fill( 0, nel, 3 )

      for il in 1:nel
          x   = readline(fid)
          xx  = split(x)

          belm[jjb][il,1] = parse(Int64,xx[2]) + 1
          belm[jjb][il,2] = parse(Int64,xx[3]) + 1
          belm[jjb][il,3] = jjb
      end
  end

  close(fid)

  bel = fill( 0, sum(nelm), 3 )
  cnt = 1
  for jjb in 1:nm
      bel[cnt:(cnt+nelm[jjb]-1),:] = belm[jjb]
      cnt += nelm[jjb]
  end

  return (p, t, bel, sm)

end

"""
    readSU2_3D( flname::String )

Reads in 3D SU2 mesh.
"""
function readSU2_3D( flname::String )

    dim = 3

    # Open file
    fid = open( flname, "r" )

    # Check dimension
    x = "%"
    while x[1] == '%'
        x = readline(fid)
    end
    xx   = split(x)
    ndim = parse( Int64, xx[end] )
    if ndim != dim
        error(" Input mesh error: Input mesh is wrong dimension, should be 3D")
    end

    # Read number of elements and triangle connectivity
    x = "%"
    while x[1] == '%'
        x = readline(fid)
    end
    xx    = split( x )
    nelem = parse( Int64, xx[end] )
    t = fill( 0, nelem, 4 )
    for il in 1:nelem
        x   = readline(fid)
        xx  = split(x)
        ind = parse(Int64,xx[end]) + 1

        t[ind,1] = parse(Int64,xx[2]) + 1
        t[ind,2] = parse(Int64,xx[3]) + 1
        t[ind,3] = parse(Int64,xx[4]) + 1
        t[ind,4] = parse(Int64,xx[5]) + 1
    end
    if parse(Int64,xx[1]) != 10
        error(" Input mesh error: Only tetrahedra are supported for 3D")
    end

    # Read node locations
    x = "%"
    while x[1] == '%'
        x = readline(fid)
    end
    xx = split( x )
    np = parse( Int64, xx[end] )
    p  = fill( 0.0, np, 3 )
    for il in 1:np
        x   = readline(fid)
        xx  = split(x)
        ind = parse(Int64,xx[end]) + 1

        p[ind,1] = parse(Float64,xx[1])
        p[ind,2] = parse(Float64,xx[2])
        p[ind,3] = parse(Float64,xx[3])
    end

    # Read boundary elements
    x = "%"
    while x[1] == '%'
        x = readline(fid)
    end
    xx = split( x )
    nm = parse( Int64, xx[end] )
    nnm = fill( 0, nm )
    sm = fill( "tag", nm )
    belm = fill( fill(0, 0, 0), nm )
    nelm = fill( 0, nm )

    for jjb in 1:nm
        x         = readline( fid )
        xx        = split( x )
        sm[jjb]   = xx[end]

        x         = readline( fid )
        xx        = split( x )
        nel       = parse( Int64, xx[end] )
        nelm[jjb] = nel
        belm[jjb] = fill( 0, nel, 4 )

        for il in 1:nel
            x   = readline(fid)
            xx  = split(x)

            belm[jjb][il,1] = parse(Int64,xx[2]) + 1
            belm[jjb][il,2] = parse(Int64,xx[3]) + 1
            belm[jjb][il,3] = parse(Int64,xx[4]) + 1
            belm[jjb][il,4] = jjb
        end
    end

    close(fid)

    bel = fill( 0, sum(nelm), 4 )
    cnt = 1
    for jjb in 1:nm
        bel[cnt:(cnt+nelm[jjb]-1),:] = belm[jjb]
        cnt += nelm[jjb]
    end

    return (p, t, bel, sm)

end
