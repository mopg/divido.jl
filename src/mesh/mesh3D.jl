type Mesh3D <: Mesh

  porder::Int64     # Polynomial order

  n::Int64          # Number of nodes

  p::Array{Float64} # Nodal locations
  t::Array{Int64}   # Triangle - node connectivity
  t2f::Array{Int64} # Triangle - face connectivity
  f::Array{Int64}   # Face - node/triangle connectivity

end
