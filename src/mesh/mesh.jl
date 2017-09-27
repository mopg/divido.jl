# ---------------------------------------------------------------------------- #
#
#   mesh.jl
#
#   Abstract mesh type
#   This allows for writing an n-dimensional version of solver methods
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

abstract type Mesh

end

using GroupSlices # TODO: KEEP TRACK OF DEPENDENCIES

include("mesh2D.jl")
include("mesh3D.jl")
