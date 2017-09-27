using GroupSlices # TODO: KEEP TRACK OF DEPENDENCIES

abstract type Mesh

end

include("mesh2D.jl")
include("mesh3D.jl")

getdim(mesh::Mesh) = mesh.porder
