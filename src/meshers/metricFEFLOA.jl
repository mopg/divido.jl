"""
    writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Vector{Float64}, flname::String )

Writes 2D isotropic metric to feflo.a format in `flname`.
"""
function writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Vector{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != length(metric)
        error(" Isotropic metric error: size of metric does not correspond with size of mesh")
    end

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf( fid, "Dimension 2\n\n" )

    @printf( fid, "SolAtVertices\n" )

    @printf( fid, "%i\n", np )
    @printf( fid, "1 1\n" )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\n", metric[jj])
    end

    close(fid)

    # TODO: Corners

end

"""
    writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Matrix{Float64}, flname::String )

Writes 2D anisotropic metric to feflo.a format in `flname`.
"""
function writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Matrix{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != size(metric,1)
        error(" Ansotropic metric error: size of metric does not correspond with size of mesh")
    end
    if size(metric,2) != 3
        error(" Ansotropic metric error: wrong number of entries for metric")
    end

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf( fid, "Dimension 2\n\n" )

    @printf( fid, "SolAtVertices\n" )

    @printf( fid, "%i\n", np )
    @printf( fid, "1 1\n" )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\t%16.15e\t%16.15e\n", metric[jj,1], metric[jj,2], metric[jj,3])
    end

    close(fid)

    # TODO: Corners

end

"""
    writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Vector{Float64}, flname::String )

Writes 3D isotropic metric to feflo.a format in `flname`.
"""
function writeMetric( mesh::Mesh3D, mesher::MesherFEFLOA,
                      metric::Vector{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != length(metric)
        error(" Isotropic metric error: size of metric does not correspond with size of mesh")
    end

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf( fid, "Dimension 3\n\n" )

    @printf( fid, "SolAtVertices\n" )

    @printf( fid, "%i\n", np )
    @printf( fid, "1 1\n" )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\n", metric[jj])
    end

    close(fid)

    # TODO: Corners

end

"""
    writeMetric( mesh::Mesh2D, mesher::MesherFEFLOA,
                      metric::Matrix{Float64}, flname::String )

Writes 3D anisotropic metric to feflo.a format in `flname`.
"""
function writeMetric( mesh::Mesh3D, mesher::MesherFEFLOA,
                      metric::Matrix{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != size(metric,1)
        error(" Ansotropic metric error: size of metric does not correspond with size of mesh")
    end
    if size(metric,2) != 3
        error(" Ansotropic metric error: wrong number of entries for metric")
    end

    fid = open( flname, "w" )

    @printf( fid, "MeshVersionFormatted 2\n\n" )

    @printf( fid, "Dimension 3\n\n" )

    @printf( fid, "SolAtVertices\n" )

    @printf( fid, "%i\n", np )
    @printf( fid, "1 1\n" )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\t%16.15e\t%16.15e\t%16.15e\t%16.15e\t%16.15e\n",
         metric[jj,1], metric[jj,2], metric[jj,3],
         metric[jj,4], metric[jj,5], metric[jj,6])
         # Ordered (1,1), (2,1), (2,2), (3,1), (3,2), (3,3)
    end

    close(fid)

    # TODO: Corners

end
