# ---------------------------------------------------------------------------- #
#
#   mesher.jl
#
#   Abstract mesher type
#   This allows for writing an n-dimensional version of mesher methods
#
#   λυτέος
#   Fall 2017
#
#   Max Opgenoord
#
# ---------------------------------------------------------------------------- #

"""
    Mesher

Mesher abstract type:
Overarching abstract type for meshers (BAMG, FEFLOA).
"""
abstract type Mesher

end

include("mesherBAMG.jl")
include("mesherFEFLOA.jl")

include("metricBAMG.jl")
include("metricFEFLOA.jl")
