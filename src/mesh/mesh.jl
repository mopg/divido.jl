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

"""
    Mesh

Mesh abstract type:
Overarching abstract type for mesh types (2D and 3D).
"""
abstract type Mesh

end

using GroupSlices # TODO: KEEP TRACK OF DEPENDENCIES

include("mesh2D.jl")
include("mesh3D.jl")
