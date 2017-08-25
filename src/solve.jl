type Solve

  p::Int64 # Number of states

  mesh::Mesh # dynamic constraints
  mat::Material

end

function Solve( nx )

  model( nx, emptyfunc, emptyfunc, emptyfunc, emptyfunc, emptyfunc )

end

function print( )

  println( "bla" )

end
