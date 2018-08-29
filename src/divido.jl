# ---------------------------------------------------------------------------- #
#
#   divido.jl
#
#   Mesh data structures for FEM and lattice structures
#
#   divido
#   Fall 2017 / Spring 2018
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

__precompile__()

"""
    divido

Julia package for mesh data structures for FEM and lattice structures.

Max Opgenoord

Fall 2017 / Spring 2018
"""

module divido

# mesh
export Mesh2D, Mesh3D
include("mesh/mesh.jl")
include("mesh/master.jl")
include("mesh/compJacob.jl")

# links to meshers
export writeMsh, readBAMG, runMesher, MesherFEFLOA, MesherBAMG, writeMsh
export readFEFLOA_2D, readFEFLOA_3D, writeMetric
include("meshers/mesher.jl")

# io
include("io/readSU2.jl")

end
