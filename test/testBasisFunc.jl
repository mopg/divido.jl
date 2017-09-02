
# order 1
order = 2
x = [1.0, 2.0]
y = [2.0, 1.0]
(phi, dphi) = basisFuncTriangleLeg( Val{order}, x, y )

@test phi[1,1] == 1.0

@test abs(dphi[1,1,1]) < 1e-15
@test abs(dphi[1,1,2]) < 1e-15
