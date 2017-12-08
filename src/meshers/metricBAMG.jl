"""
    writeMetric( mesh::Mesh2D, mesher::MesherBAMG,
                      metric::Vector{Float64}, flname::String )

Writes 2D isotropic metric to BAMG format in `flname`.
"""
function writeMetric( mesh::Mesh2D, mesher::MesherBAMG,
                      metric::Vector{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != length(metric)
        error(" Isotropic metric error: size of metric does not correspond with size of mesh")
    end

    fid = open( flname, "w" )

    @printf( fid, "%i 1\n", np )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\n", metric[jj])
    end

    close(fid)

    # TODO: Corners

end

"""
    writeMetric( mesh::Mesh2D, mesher::MesherBAMG,
                      metric::Matrix{Float64}, flname::String )

Writes 2D anisotropic metric to BAMG format in `flname`.
"""
function writeMetric( mesh::Mesh2D, mesher::MesherBAMG,
                      metric::Matrix{Float64}, flname::String )

    np = size(mesh.p,1)
    if np != size(metric,1)
        error(" Ansotropic metric error: size of metric does not correspond with size of mesh")
    end
    if size(metric,2) != 3
        error(" Ansotropic metric error: wrong number of entries for metric")
    end

    fid = open( flname, "w" )

    @printf( fid, "%i 3\n", np )

    # Vertices
    for jj = 1:np
        @printf(fid, "%16.15e\t%16.15e\t%16.15e\n", metric[jj,1], metric[jj,2], metric[jj,3])
    end

    close(fid)

    # TODO: Corners

end
