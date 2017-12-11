# ---------------------------------------------------------------------------- #
#
#   mesherBAMG.jl
#
#   BAMG-specific read and write
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    MesherBAMG

MesherBAMG type:
Type for BAMG meshes.
"""
struct MesherBAMG <: Mesher

end

"""
    writeMsh( mesh::Mesh2D, mesher::MesherBAMG, flname::String )

Writes 2D mesh in `mesh` to `flname` in format used by BAMG.
"""
function writeMsh( mesh::Mesh2D, mesher::MesherBAMG, flname::String )

    np = size(mesh.p,1)

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 0\n\n" )

    @printf(fid, "Dimension 2\n\n")

    # Vertices
    @printf(fid, "Vertices\n%i\n", np)
    for jj = 1:size(mesh.p,1)
        @printf(fid, "%16.15e\t%16.15e\t%i\n", mesh.p[jj,1], mesh.p[jj,2], 1)
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

    # Triangles
    @printf(fid, "Triangles\n%i\n",size(mesh.t,1))
    for jj = 1:size(mesh.t,1)
        @printf(fid, "%i\t%i\t%i\t%i\n", mesh.t[jj,1], mesh.t[jj,2], mesh.t[jj,3], 1)
    end
    @printf(fid,"\n")

    close(fid)

    # TODO: Corners

end

"""
    readBAMG( flname::String )

Reads (2D) BAMG mesh from `flname`.
"""
function readBAMG( flname::String )

    dim = 2

    # Open file
    fid = open( flname, "r" )

    # Check dimension
    for il in 1:3
        x = readline(fid)
    end
    x = readline(fid)
    if parse(Int64, x) != dim
        error(" BAMG read-in: dimension should be 2")
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

    close(fid)

    return( p, t, bel )

end
