# top-k-sum
Routines for Euclidean projection of an $n$-dimensional vector $x^0$ onto the top- $\\!\\!k$-sum constraint
```math
\mathcal{B}_{(k)}^r \coloneqq \Bigl\{x\in\mathbb{R}^n : \sum_{i=1}^k x_{[i]} \leq r\Bigr\},
```
where $k\in\\{1,\ldots,n\\}$ and $r\in\mathbb{R}$ are parameters and where $x_{[i]}$ denotes the $i^{\text{th}}$ largest element of $x$.

The repository contains two directories:
- `src/`: implementations of the algorithms summarized in XYZ (arxiv).
- `run/`: scripts for running the experiments summarized in XYZ (arxiv). Experiments can be run from the `run/` folder in terminal as `$ julia run_tks_experiments.jl` by specifying a proper `DATAPATH` variable.

Routines in `src/` include three finite-termination algorithms:
- ESGS: an early-stopping grid-search approach callable by `project_topksum_esgs!`
- PLCP: a parametric-LCP approach callable by `project_topksum_plcp!`
- GRID: a full grid-search approach callable by `project_topksum_grid!`

and two inexact, QP solvers:
- GRBS: a QP barrier method using Gurobi based on the _sorted_ formulation, and callable by `project_topksum_grbs`
- GRBU: a QP barrier method using Gurobi based on the _unsorted_ formulation, and callable by `project_topksum_grbu`. **Note**: there is not a native `quickselect` procedure in `Julia`, so GRBU assumes that the largest $k$ elements of the input vector are partially sorted

An example calling sequence for ESGS and the Gurobi-based methods is provided below.

```
using Random
Random.seed!(1234567)

n = 10;
k = 4;
r = 0.1234567;
x0 = randn(n);
pi = sortperm(x0, rev=true); # full sort
pi_partial = sortperm(x0, alg=PartialQuickSort(1:k), rev=true); # partial sort
x0sort = x0[pi]; # full sort
x0psort = x0[pi_partial]; # partial sort

active = sum(x0sort[1:k]) > r; # binary: is the constraint active
xbarsort = similar(x0sort);
xbarsort .= x0sort;
prepop = true; # binary: is the solution vector pre-populated with elements from sorted input vector
out_esgs = project_topksum_esgs!(xbarsort, x0sort, r, k, active, prepop); # solve with ESGS method
out_grbs = project_topksum_grbs(x0sort, r, k); # solve with GRBS method
out_grbu = project_topksum_grbu(x0psort, r, k); # solve with GRBU (partial sort) method
```
which results in
```
julia> # display sorted solution
       vcat(["sorted input" "esgs" "grbs" "grbu"], hcat(x0sort, xbarsort, out_grbs[:x], out_grbu[:xsort]))
11×4 Matrix{Any}:
   "sorted input"    "esgs"     "grbs"     "grbu"
  2.0882            0.762312   0.762312   0.762312
  1.01192          -0.212952  -0.212952  -0.212952
  1.00237          -0.212952  -0.212952  -0.212952
  0.577227         -0.212952  -0.212952  -0.212952
  0.377373         -0.212952  -0.212952  -0.212952
 -0.0814697        -0.212952  -0.212952  -0.212952
 -0.187464         -0.212952  -0.212952  -0.212952
 -0.47593          -0.47593   -0.47593   -0.47593
 -0.754928         -0.754928  -0.754928  -0.754928
 -1.41154          -1.41154   -1.41154   -1.41154

julia> # verify feasibility
       (sum(xbarsort[1:k]), sum(out_grbs[:x][1:k]), sum(out_grbu[:xsort][1:k]), r)
(0.12345669999999945, 0.12345669999491526, 0.12345669999960365, 0.1234567)

julia> # display objective value
       objvals = (0.5sum((xbarsort .- x0sort).^2), 0.5sum((out_grbs[:x] .- x0sort).^2), 0.5sum((out_grbu[:xsort] .- x0sort).^2));
julia> minobj = minimum(objvals);
julia> objvals .- minobj
(0.0, 2.2381181352670865e-9, 2.375699637013895e-10)

julia> # display index-pairs (k0,k1)
       (out_esgs[2], get_k0k1(out_grbs[:x], k), get_k0k1(out_grbu[:xsort], k))
((1, 7), (3, 4), (3, 4))
```
Note that the Gurobi-based methods return incorrect index-pairs, even in an $n=10$ dimensional problem.
In general, it is not obvious how to recover appropriate $(k_0,k_1)$ information from a Gurobi solution.
