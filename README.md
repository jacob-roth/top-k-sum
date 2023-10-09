# top-k-sum
Routines for Euclidean projection of an $n$-dimensional vector $x^0$ onto the top- $\\!\\!k$-sum constraint
```math
\mathcal{B}_{(k)}^r \coloneqq \Bigl\{x\in\mathbb{R}^n : \sum_{i=1}^k x_{[i]} \leq r\Bigr\},
```
where $k\in\\{1,\ldots,n\\}$ and $r\in\mathbb{R}$ are parameters and where $x_{[i]}$ denotes the $i^{\text{th}}$ largest element of $x$.

The repository contains two directories:
- `src/`: implementations of the algorithms summarized in XYZ (arxiv).
- `run/`: scripts for running the experiments summarized in XYZ (arxiv). Experiments can be run from the `run/` folder in terminal as `$ julia run_tks_experiments.jl`.

Routines in `src/` include three finite-termination algorithms:
- ESGS: an early-stopping grid-search approach callable by `project_topksum_esgs!`
- PLCP: a parametric-LCP approaach callable by `project_topksum_plcp!`
- GRID: a full grid-search approach callable by `project_topksum_grid!`

and two inexact, QP solvers:
- GRBS: a QP barrier method using Gurobi based on the _sorted_ formulation, and callable by `project_topksum_gurobi`
- GRBU: a QP barrier method using Gurobi based on the _unsorted_ formulation, and callable by `project_topksum_unsort_gurobi`

An example calling sequence for ESGS is provided below.

```
using Random
Random.seed!(1234567)

n = 10;
k = 4;
r = 0.1234567;
x0 = randn(n);
pi = sortperm(x0, rev=true);
pi_inv = invperm(pi);
x0sort = x0[pi];

xbarsort = similar(x0sort);
xbarsort .= x0sort;
active = sum(x0sort[1:k]) > r;
prepop = true; # solution vector is pre-populated with elements from sorted input vector
out = project_maxksum_esgs!(xbarsort, x0sort, r, k, active, prepop); # solve with ESGS method
```
which results in
```
julia> # display sorted solution
       vcat(["sorted input" "sorted solution"], hcat(x0sort, xbarsort))
11×2 Matrix{Any}:
   "sorted input"    "sorted solution"
  2.0882            0.762312
  1.01192          -0.212952
  1.00237          -0.212952
  0.577227         -0.212952
  0.377373         -0.212952
 -0.0814697        -0.212952
 -0.187464         -0.212952
 -0.47593          -0.47593
 -0.754928         -0.754928
 -1.41154          -1.41154

julia> # display original solution
       vcat(["unsorted input" "unsorted solution"], hcat(x0, xbarsort[pi_inv]))
11×2 Matrix{Any}:
   "unsorted input"    "unsorted solution"
  0.577227           -0.212952
 -0.187464           -0.212952
  1.00237            -0.212952
  1.01192            -0.212952
  2.0882              0.762312
  0.377373           -0.212952
 -0.0814697          -0.212952
 -1.41154            -1.41154
 -0.47593            -0.47593
 -0.754928           -0.754928

julia> # verify feasibility
       (sum(xbarsort[1:k]), r)
(0.12345669999999945, 0.1234567)
```
