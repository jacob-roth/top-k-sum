using CxxWrap
module CppOWL
  using CxxWrap
  import CxxWrap: StdVector

  # Define a function that returns the path to the dynamic library
  function get_lib_path()
    # Ensure this exactly matches the library name and path
    return joinpath(@__DIR__, "CppOWL/build/libowlball.dylib")  # Correct path and extension
  end

  # Correctly pass the function to @wrapmodule
  @wrapmodule(get_lib_path)  # Note the absence of ()

  # # Julia convenience wrapper for evaluateProx that handles vector conversions
  # function evaluateProx(z_in::AbstractVector{Tf}, w::AbstractVector{Tf}, epsilon::Tf, x_out::AbstractVector{Tf}, sorted_and_positive::Bool) where {Tf<:AbstractFloat}
  #   # Call the C++ function
  #   evaluateProx(StdVector{Tf}(z_in), StdVector{Tf}(w), epsilon, StdVector{Tf}(x_out), sorted_and_positive)
  # end

  function __init__()
    @initcxx
  end
end
