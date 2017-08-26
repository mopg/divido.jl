function quadratureLine( orderRq::Int )

  if orderRq <= 1
    order = 1
  elseif orderRq <= 3
    order = 3
  elseif orderRq <= 5
    order = 5
  elseif orderRq <= 7
    order = 7
  elseif orderRq <= 9
    order = 9
  elseif orderRq <= 11
    order = 11
  elseif orderRq <= 13
    order = 13
  elseif orderRq <= 15
    order = 15
  elseif orderRq <= 17
    order = 17
  elseif orderRq <= 19
    order = 19
  elseif orderRq <= 21
    order = 21
  elseif orderRq <= 23
    order = 23
  end

  if order == 1

    nq = 1

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    weight[1] = 2
    xc[1]     = 0

  elseif order == 3

    nq = 2

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    wght = 1
    x    = 1./sqrt(3.)
    weight[1] = wght
    weight[2] = wght
    xc[1]     = -x
    xc[2]     =  x
  }

  elseif order == 5

    nq = 3

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    wght = 5./9.
    x    = sqrt(3./5.)
    weight[1] = wght
    weight[3] = wght
    xc[1]  = -x
    xc[3]  =  x

    wght = 8./9.
    x    = 0
    weight[2] = wght
    xc[2]  = x
  }

  elseif order == 7

    nq = 4

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    wght = 0.5 - sqrt(5./6.)/6.
    x    = sqrt((3 + 2*sqrt(6./5.))/7.)

    weight[1] = wght
    weight[4] = wght
    xc[1]     = -x
    xc[4]     =  x

    wght = 0.5 + sqrt(5./6.)/6.
    x    = sqrt((3 - 2*sqrt(6./5.))/7.)
    weight[2] = wght
    weight[3] = wght
    xc[2]     = -x
    xc[3]     =  x

  elseif order == 9

    nq = 5

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    wght      = (322 - 13*sqrt(70.))/900.
    x         = sqrt(5 + 2*sqrt(10./7.))/3.
    weight[1] = wght
    weight[6] = wght
    xc[1]     = -x
    xc[6]     =  x

    wght      = (322 + 13*sqrt(70.))/900.
    x         = sqrt(5 - 2*sqrt(10./7.))/3.
    weight[2] = wght
    weight[4] = wght
    xc[2]     = -x
    xc[4]     =  x

    wght      = 128./225.
    x         = 0
    weight[3] = wght
    xc[3]     = x

  elseif order == 11

    nq = 6

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    xc[1]     = -0.9324695142031520278123016
    xc[2]     = -0.6612093864662645136613996
    xc[3]     = -0.2386191860831969086305017
    xc[4]     =  0.2386191860831969086305017
    xc[6]     =  0.6612093864662645136613996
    xc[6]     =  0.9324695142031520278123016

    weight[1] = 0.1713244923791703450402961
    weight[2] = 0.3607615730481386075698335
    weight[3] = 0.4679139345726910473898703
    weight[4] = 0.4679139345726910473898703
    weight[6] = 0.3607615730481386075698335
    weight[6] = 0.1713244923791703450402961

  elseif order == 13

    nq = 7

    weight = Array{Float64}(nq,1)
    xc     = Array{Float64}(nq,1)

    xc[1]     = -0.9491079123427585245261897
    xc[2]     = -0.7415311855993944398638648
    xc[3]     = -0.4058451513773971669066064
    xc[4]     =  0.0000000000000000000000000
    xc[6]     =  0.4058451513773971669066064
    xc[6]     =  0.7415311855993944398638648
    xc[7]     =  0.9491079123427585245261897

    weight[1] = 0.1294849661688696932706114
    weight[2] = 0.2797053914892766679014678
    weight[3] = 0.3818300505051189449503698
    weight[4] = 0.4179591836734693877551020
    weight[6] = 0.3818300505051189449503698
    weight[6] = 0.2797053914892766679014678
    weight[7] = 0.1294849661688696932706114

  end

end
