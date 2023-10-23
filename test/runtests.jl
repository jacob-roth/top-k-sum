# using TopKSum
include("../src/TopKSum.jl")
using Test
using LinearAlgebra
using Gurobi, JuMP
using SparseArrays, SparseMatricesCSR
using Random

#=
small problem
=#
@testset "small problem" begin
Random.seed!(1234567)
tol = 1e-7;
n = 10;
k = 4;
r = 0.1234567;
x0 = randn(n);
pi = sortperm(x0, rev=true); # full sort
pi_partial = sortperm(x0, alg=PartialQuickSort(1:k), rev=true); # partial sort
x0sort = x0[pi]; # full sort
x0psort = x0[pi_partial]; # partial sort

active = sum(x0sort[1:k]) > r; # binary: is the constraint active
xbarsort_esgs = similar(x0sort); xbarsort_esgs .= x0sort;
xbarsort_grid = similar(x0sort); xbarsort_grid .= x0sort;
xbarsort_plcp = similar(x0sort); xbarsort_plcp .= x0sort;
prepop = true; # binary: is the solution vector pre-populated with elements from sorted input vector
out_esgs = project_topksum_esgs!(xbarsort_esgs, x0sort, r, k, active, prepop); # solve with ESGS method
out_grid = project_topksum_grid!(xbarsort_grid, x0sort, r, k, active, prepop); # solve with ESGS method
out_plcp = project_topksum_plcp!(xbarsort_plcp, x0sort, r, k, active, prepop); # solve with ESGS method
out_grbs = project_topksum_grbs(x0sort, r, k); # solve with GRBS method
out_grbu = project_topksum_grbu(x0psort, r, k); # solve with GRBU (partial sort) method
@test norm(xbarsort_esgs .- out_grbs[:x], Inf) <= tol
@test norm(xbarsort_esgs .- out_grbu[:xsort], Inf) <= tol
@test norm(xbarsort_esgs .- xbarsort_grid, Inf) <= tol
@test norm(xbarsort_esgs .- xbarsort_plcp, Inf) <= tol
@test sum(xbarsort_esgs[1:k]) <= r+eps()
@test sum(out_grbs[:x][1:k]) <= r+eps()
@test sum(out_grbu[:xsort][1:k]) <= r+eps()
end # small problem

#=
constraint matrices
=#
function project_polyhedron(x0::AbstractVector, B::AbstractMatrix, b::AbstractVector)
  """project onto {x : Bx<=b}"""
  n = length(x0)
  model = Model(Gurobi.Optimizer)
  @variable(model, x[1:n])
  @constraint(model, B*x .<= b)
  @objective(model, Min, 0.5 * sum( (x0 .- x).^2 ))
  optimize!(model)
  return value.(x)
end
@testset "constraint matrices" begin
Random.seed!(1234567)
tol = 1e-6;
n = 10;
k = 4;
r = 0.1234567;
x0 = randn(n);
pi = sortperm(x0, rev=true); # full sort
pi_partial = sortperm(x0, alg=PartialQuickSort(1:k), rev=true); # partial sort
x0sort = x0[pi]; # full sort
x0psort = x0[pi_partial]; # partial sort

active = sum(x0sort[1:k]) > r; # binary: is the constraint active
xbarsort = similar(x0sort); xbarsort .= x0sort;
prepop = true; # binary: is the solution vector pre-populated with elements from sorted input vector
out_esgs = project_topksum_esgs!(xbarsort, x0sort, r, k, active, prepop); # solve with ESGS method

#
# primal
#

# isotonic: eq. 7
B = get_isotonic(n, k)
b = [r; zeros(n-1)]
xtilde = project_polyhedron(x0sort, B, b)
xbar = sort(xtilde, rev=true)
@test norm(xbar-xbarsort, Inf) <= tol

# unsorted-top-k: eq. 5
B = get_unsortedtopk(x0, k)
b = [r; zeros(n-1)]
xtilde = project_polyhedron(x0, B, b)
xbar = sort(xtilde, rev=true)
@test norm(xbar-xbarsort, Inf) <= tol

# partial-sorted-top-k: (based on eq. 5)
B = get_unsortedtopk(x0psort, k)
b = [r; zeros(n-1)]
xtilde = project_polyhedron(x0psort, B, b)
xbar = sort(xtilde, rev=true)
@test norm(xbar-xbarsort, Inf) <= tol

#
# dual (moreau)
#

# isotonic: eq. 7
Binvt = get_isotonic_moreau(n, k) # == inv(get_isotonic(n, k))'
ind = ones(n);
ybar = project_polyhedron(x0sort .- r/k .* ind, -Binvt, zeros(n))
xbar = x0sort .- ybar
@test norm(xbar-xbarsort, Inf) <= 2tol
@test norm(inv(get_isotonic(n, k))' - Binvt) == 0

# unsorted-top-k: eq. 6
Binvt = Matrix(inv(get_unsortedtopk(x0, k))')
ind = ones(n);
ybar = project_polyhedron(x0 .- r/k .* ind, -Binvt, zeros(n))
xbar = sort(x0 .- ybar, rev=true)
@test norm(xbar-xbarsort, Inf) <= tol

# partial-sorted-top-k: eq. 6
Binvt = Matrix(inv(get_unsortedtopk(x0psort, k))')
ind = ones(n);
ybar = project_polyhedron(x0sort .- r/k .* ind, -Binvt, zeros(n))
xbar = x0psort .- ybar
@test norm(xbar-xbarsort, Inf) <= 2tol

end # constraint matrices


#===========commented-out===========
#=
lcp matrices
=#
@testset "lcp matrices" begin
Random.seed!(1234567)
tol = 1e-6;
n = 10;
k = 4;
r = 0.1234567;
x0 = randn(n);
pi = sortperm(x0, rev=true); # full sort
pi_partial = sortperm(x0, alg=PartialQuickSort(1:k), rev=true); # partial sort
x0sort = x0[pi]; # full sort
x0psort = x0[pi_partial]; # partial sort

active = sum(x0sort[1:k]) > r; # binary: is the constraint active
xbarsort = similar(x0sort); xbarsort .= x0sort;
prepop = true; # binary: is the solution vector pre-populated with elements from sorted input vector
out_esgs = project_topksum_esgs!(xbarsort, x0sort, r, k, active, prepop); # solve with ESGS method

nits = [project_topksum_esgs!(xbarsort, x0sort, r, k, sum(x0sort[1:k]) > r, false)[3] for k in 1:100]
nits = [project_topksum_eric!(xbarsort, x0sort, r, k)[2] for k in 1:100]

#
# primal
#

# isotonic: eq. 7
B = get_isotonic(n, k)
M = B * B'
@test all([i==j ? true : M[i,j] <= 0 for i in 1:n for j in 1:n])
q = [r; zeros(n-1)]
Bplcp = get_isotonic(n, k)[2:end,:]
Mplcp = Bplcp * Bplcp'
qplcp = [r; zeros(n-2)]
@test all([i==j ? true : Mpclp[i,j] <= 0 for i in 1:n-1 for j in 1:n-1])

# unsorted-top-k: eq. 5
B = get_unsortedtopk(x0, k)
M = B * B'
q = [r; zeros(n-1)] - B*x0
sol = solve!(LCP(M, q))
xbarsort .- sort(x0 - B' * sol.sol, rev=true)
#? idea: what if we kept pivoting on (smaller but notsorted) entries after k; would get T_{(k)}(x) < r; is this useful?

# partial-sorted-top-k: (based on eq. 5)
B = get_unsortedtopk(x0psort, k)

#
# dual
#

# isotonic lcp: eq. 7
Binvt = get_isotonic_moreau(n, k)
M = Binvt * Binvt'
end # lcp matrices

#=
heap
=#

n = 20
k = 6
x0 = randn(n);
xbar = similar(x0);
x0sort = sort(x0, rev=true)
xbarsort = similar(x0);
r = 0.1234567
project_topksum_esgs_heap!(xbar, x0, r, k)
project_topksum_esgs!(xbarsort, x0sort, r, k, true, false)
hcat(xbarsort, sort(xbar, rev=true))
===========commented-out===========#