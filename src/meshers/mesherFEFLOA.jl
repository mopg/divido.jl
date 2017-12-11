# ---------------------------------------------------------------------------- #
#
#   mesherFEFLOA.jl
#
#   feflo.a-specific read and write
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    MesherFEFLOA

MesherFEFLOA type:
Type for FEFLOA meshes.
"""
struct MesherFEFLOA <: Mesher

end

"""
    writeMsh( mesh::Mesh2D, mesher::MesherFEFLOA, flname::String )

Writes 2D mesh in `mesh` to `flname` in format used by feflo.a.
"""
function writeMsh( mesh::Mesh2D, mesher::MesherFEFLOA, flname::String )

    np = size(mesh.p,1)

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf(fid, "Dimension 2\n\n")

    # Vertices
    @printf(fid, "Vertices\n%i\n", np)
    for jj = 1:size(mesh.p,1)
        @printf(fid, "%16.15e\t%16.15e\t%i\n", mesh.p[jj,1], mesh.p[jj,2], 1)
    end
    @printf(fid,"\n")

    # Triangles
    @printf(fid, "Triangles\n%i\n",size(mesh.t,1))
    for jj = 1:size(mesh.t,1)
        @printf(fid, "%i\t%i\t%i\t%i\n", mesh.t[jj,1], mesh.t[jj,2], mesh.t[jj,3], 1)
    end
    @printf(fid,"\n")

    # Edges
    nbound = size(mesh.fb,1)
    @printf(fid, "Edges\n%i\n",nbound)
    for ii = 1:nbound
        indf = mesh.fb[ii,1]
        @printf(fid, "%i\t%i\t%i\n", mesh.f[indf,1], mesh.f[indf,2], mesh.fb[ii,2])
    end
    @printf(fid,"\n")

    close(fid)

    # TODO: Corners

end

"""
    writeMsh( mesh::Mesh3D, mesher::MesherFEFLOA, flname::String )

Writes 3D mesh in `mesh` to `flname` in format used by feflo.a.
"""
function writeMsh( mesh::Mesh3D, mesher::MesherFEFLOA, flname::String )

    np = size(mesh.p,1)

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf(fid, "Dimension 3\n\n")

    # Vertices
    @printf(fid, "Vertices\n%i\n", np)
    for jj = 1:size(mesh.p,1)
        @printf(fid, "%16.15e\t%16.15e\t%16.15e\t%i\n", mesh.p[jj,1], mesh.p[jj,2],mesh.p[jj,3], 1)
    end
    @printf(fid,"\n")

    # Triangles
    @printf(fid, "Tetrahedra\n%i\n",size(mesh.t,1))
    for jj = 1:size(mesh.t,1)
        @printf(fid, "%i\t%i\t%i\t%i\t%i\n", mesh.t[jj,1], mesh.t[jj,2], mesh.t[jj,3], mesh.t[jj,4], 1)
    end
    @printf(fid,"\n")

    # Edges
    nbound = size(mesh.fb,1)
    @printf(fid, "Triangles\n%i\n",nbound)
    for ii = 1:nbound
        indf = mesh.fb[ii,1]
        @printf(fid, "%i\t%i\t%i\t%i\n", mesh.f[indf,1], mesh.f[indf,2], mesh.f[indf,3], mesh.fb[ii,2])
    end
    @printf(fid,"\n")

    close(fid)

    # TODO: Corners

end

"""
    readFEFLOA_2D( flname::String )

Reads 2D feflo.a mesh from `flname`.
"""
function readFEFLOA_2D( flname::String )

    dim = 2

    # Open file
    fid = open( flname, "r" )

    # Check dimension
    x = readline(fid)
    while isempty(x) || x[1] != 'D'
        x = readline(fid)
    end
    xx = split( x )
    if parse(Int64, xx[end]) != dim
        error(" FEFLOA read-in: dimension should be 2")
    end

    # Read vertices
    while isempty(x) || x[1] != 'V'
        x = readline(fid)
    end
    x  = readline(fid)
    np = parse( Int64, x )
    p  = fill( 0.0, np, 2 )
    for il in 1:np
        x  = readline(fid)
        xx = split(x)
        p[il,1] = parse( Float64, xx[1] )
        p[il,2] = parse( Float64, xx[2] )
    end

    # Read triangles
    while isempty(x) || x[1] != 'T'
        x = readline(fid)
    end
    x  = readline(fid)
    nt = parse( Int64, x )
    t  = fill( 0, nt, 3 )
    for il in 1:nt
        x  = readline(fid)
        xx = split(x)
        t[il,1] = parse( Int64, xx[1] )
        t[il,2] = parse( Int64, xx[2] )
        t[il,3] = parse( Int64, xx[3] )
    end

    # Read boundary
    while isempty(x) || x[1] != 'E'
        x = readline(fid)
    end
    x   = readline(fid)
    nb  = parse( Int64, x )
    bel = fill( 0, nb, 3 )
    for il in 1:nb
        x  = readline(fid)
        xx = split(x)
        bel[il,1] = parse( Int64, xx[1] )
        bel[il,2] = parse( Int64, xx[2] )
        bel[il,3] = parse( Int64, xx[3] )
    end

    close(fid)

    return( p, t, bel )

end

"""
    readFEFLOA_3D( flname::String )

Reads 3D feflo.a mesh from `flname`.
"""
function readFEFLOA_3D( flname::String )

    dim = 3

    # Open file
    fid = open( flname, "r" )

    # Check dimension
    x = readline(fid)
    while isempty(x) || x[1] != 'D'
        x = readline(fid)
    end
    xx = split( x )
    if parse(Int64, xx[end]) != dim
        error(" FEFLOA read-in: dimension should be 3")
    end

    # Read vertices
    while isempty(x) || x[1] != 'V'
        x = readline(fid)
    end
    x  = readline(fid)
    np = parse( Int64, x )
    p  = fill( 0.0, np, 3 )
    for il in 1:np
        x  = readline(fid)
        xx = split(x)
        p[il,1] = parse( Float64, xx[1] )
        p[il,2] = parse( Float64, xx[2] )
        p[il,3] = parse( Float64, xx[3] )
    end

    # Read triangles
    while isempty(x) || x[1] != 'T'
        x = readline(fid)
    end
    x  = readline(fid)
    nt = parse( Int64, x )
    t  = fill( 0, nt, 4 )
    for il in 1:nt
        x  = readline(fid)
        xx = split(x)
        t[il,1] = parse( Int64, xx[1] )
        t[il,2] = parse( Int64, xx[2] )
        t[il,3] = parse( Int64, xx[3] )
        t[il,4] = parse( Int64, xx[4] )
    end

    # Read boundary
    while isempty(x) || x[1] != 'T'
        x = readline(fid)
    end
    x   = readline(fid)
    nb  = parse( Int64, x )
    bel = fill( 0, nb, 4 )
    for il in 1:nb
        x  = readline(fid)
        xx = split(x)
        bel[il,1] = parse( Int64, xx[1] )
        bel[il,2] = parse( Int64, xx[2] )
        bel[il,3] = parse( Int64, xx[3] )
        bel[il,4] = parse( Int64, xx[4] )
    end

    close(fid)

    return( p, t, bel )

end
