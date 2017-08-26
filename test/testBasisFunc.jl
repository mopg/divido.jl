include("../src/basisFuncTriangle.jl")

# order 1
order = 1
x = [1.0, 2.0]
y = [2.0, 1.0]
(phi, dphi) = basisFuncTriangle( order, x, y )

@test phi[1,1] == 1.0

@test dphi[1,1,1] == 0.0
@test dphi[1,1,2] == 0.0
