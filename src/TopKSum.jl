module TopKSum

using LinearAlgebra

# todo: in-place binary heap
# e.g., `BinaryHeap!{Float64,DataStructures.FasterReverse}(htree, x0)`
# include("projection_esgs_heap.jl")

include("projection_esgs.jl")
include("projection_plcp.jl")
include("projection_grid.jl")
include("projection_ssn.jl") # ssn solver and helpers prox pava
include("matrices.jl")
include("utility.jl")

using Requires
@require Gurobi begin
  include("projection_gurobi.jl")
  export project_topksum_grbs
  export project_topksum_grbu
end
# @require CxxWrap begin
#   include("cppowl.jl")
# end

export project_topksum_esgs!
export project_topksum_plcp!
export project_topksum_grid!
export project_topksum_ssn!

end # module
