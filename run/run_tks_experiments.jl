import Pkg; Pkg.activate("."); Pkg.instantiate()
using Distributed
addprocs(2)
@everywhere begin
  import Pkg; Pkg.activate("."); Pkg.instantiate()
  using LinearAlgebra, Random, DelimitedFiles
  using Gurobi, JuMP, SparseMatricesCSR, SparseArrays
  global const PROJPATH = match(r".*top-k-sum/", @__DIR__).match
  include(PROJPATH * "src/projection_esgs.jl")
  include(PROJPATH * "src/projection_plcp.jl")
  include(PROJPATH * "src/projection_grid.jl")
  include(PROJPATH * "src/projection_gurobi.jl")
  include(PROJPATH * "src/utility.jl")
  include(PROJPATH * "run/helper_tks_experiments.jl")
  include(PROJPATH * "run/projection.jl")
  # global const DATAPATH = "/home/roth0674/drive/tks_results/"
  global const DATAPATH = "/Users/jakeroth/tks_results_test/"
  using NaNStatistics
end # @everywhere

#* 2 procs at 10^7 takes ~15gb ram for finite-termination algs and 40gb ram for GRB*#

#
# experiments
#

# warmup
rlevel = [[-8; -4; -2; -1; -1//2; -1//10; 0] .// 1; [1//10, 5//10, 9//10, 99//100, 999//1000]]
klevel = [1//10_000, 1//1000, 1//100, 5//100, 1//10, 5//10, 9//10, 99//100, 999//1000, 9999//10_000]
DTs = [Float64]
nlevel = 10 .^ collect(1:1)
nrep = 100
maxn_gurobi = 10^5
maxn_grid = 10^5
Random.seed!(12345)
expers = nothing
rep_offset = 0
for ni in eachindex(nlevel)
  for DTi in eachindex(DTs)
    DT = DTs[DTi]
    n = nlevel[ni]
    println("===================")
    println("   n=$n, DT=$DT")
    println("===================")
    test_tks_timing(n, DT, nrep, rlevel, klevel, DATAPATH, maxn_gurobi, maxn_grid, expers, rep_offset)
  end
end

# most cases
rlevel = [[-8; -4; -2; -1; -1//2; -1//10; 0] .// 1; [1//10, 5//10, 9//10, 99//100, 999//1000]]
klevel = [1//10_000, 1//1000, 1//100, 5//100, 1//10, 5//10, 9//10, 99//100, 999//1000, 9999//10_000]
DTs = [Float64]
nlevel = 10 .^ collect(1:7)
nrep = 100
maxn_gurobi = 10^5
maxn_grid = 10^5
Random.seed!(12345)
expers = nothing
rep_offset = 0
for ni in eachindex(nlevel)
  for DTi in eachindex(DTs)
    DT = DTs[DTi]
    n = nlevel[ni]
    println("===================")
    println("   n=$n, DT=$DT")
    println("===================")
    test_tks_timing(n, DT, nrep, rlevel, klevel, DATAPATH, maxn_gurobi, maxn_grid, expers, rep_offset)
  end
end

#
# sort time
#

nlevel = 10 .^ collect(1:7)

function sort_time(
  nlevelj::Vector, nrep::Integer
)
  t = zeros(length(nlevelj), nrep)
  for j in eachindex(nlevelj)
    println("n = $(nlevelj[j])")
    rep_t = time()
    for rep in 1:nrep
      if mod(rep,10) == 0
        println("$(rep/nrep), s = $(round(time()-rep_t,digits=2))")
      end
      x0 = rand(nlevelj[j])
      t[j,rep] = @elapsed sort!(x0, rev=true)
    end
  end
  return t
end
st = sort_time(nlevel, 100)
writedlm(DATAPATH * "sort_time.csv", hcat(nlevel, st))

function partial_sort_time(
  nlevelj::Vector, nrep::Integer, kpct::Real=0.01
)
  t = zeros(length(nlevelj), nrep)
  for j in eachindex(nlevelj)
    println("n = $(nlevelj[j])")
    rep_t = time()
    for rep in 1:nrep
      if mod(rep,10) == 0
        println("$(rep/nrep), s = $(round(time()-rep_t,digits=2))")
      end
      x0 = rand(nlevelj[j])
      k = Int(ceil(length(x0) * kpct))
      # t[j,rep] = @elapsed partialsort!(x0, 1:k, rev=true)
      t[j,rep] = @elapsed sort!(x0, alg=PartialQuickSort(1:k), rev=true)
    end
  end
  return t
end
pst = partial_sort_time(nlevel, 100)
writedlm(DATAPATH * "partial_sort_time_1pct.csv", hcat(nlevel, pst))

#
# large experiments
#

# large cases
rlevel = [[-8; -4; -2; -1; -1//2; -1//10; 0] .// 1; [1//10, 5//10, 9//10, 99//100, 999//1000]]
klevel = [1//10_000, 1//1000, 1//100, 5//100, 1//10, 5//10, 9//10, 99//100, 999//1000, 9999//10_000]
expers = [
  (-4//1, 1//1000), (-1//10, 1//1000), (1//10, 1//1000), (99//100, 1//1000),
  (-4//1, 1//20),   (-1//10, 1//20),   (1//10, 1//20),   (99//100, 1//20),
]
DTs = [Float64]
nlevel = 10 .^ [1;6;7]
nrep = 2
maxn_gurobi = 10^7
maxn_grid = 10^7
Random.seed!(12345)
rep_offset = 100
for ni in eachindex(nlevel)
  for DTi in eachindex(DTs)
    DT = DTs[DTi]
    n = nlevel[ni]
    println("===================")
    println("   n=$n, DT=$DT")
    println("===================")
    test_tks_timing(n, DT, nrep, rlevel, klevel, DATAPATH, maxn_gurobi, maxn_grid, expers, rep_offset)
  end
end
