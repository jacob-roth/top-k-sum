function initialize_output(n, DT, nr, nk)
  alg_names = ["ETGS", "PLCP", "GRID", "GRBS", "GRBU"]
  out = Dict()
  out["ETGS"] = Dict()
  out["PLCP"] = Dict()
  out["GRID"] = Dict()
  out["GRBS"] = Dict()
  out["GRBU"] = Dict()
  out["ALL"] = Dict()
  # k0
  out["ETGS"][:k0] = fill(0, nr, nk)
  out["PLCP"][:k0] = fill(0, nr, nk)
  out["GRID"][:k0] = fill(0, nr, nk)
  out["GRBS"][:k0] = fill(0, nr, nk)
  out["GRBU"][:k0] = fill(0, nr, nk)
  
  # k1
  out["ETGS"][:k1] = fill(0, nr, nk)
  out["PLCP"][:k1] = fill(0, nr, nk)
  out["GRID"][:k1] = fill(0, nr, nk)
  out["GRBS"][:k1] = fill(0, nr, nk)
  out["GRBU"][:k1] = fill(0, nr, nk)
  
  # t_init
  out["ETGS"][:t_init] = fill(+Inf, nr, nk)
  out["PLCP"][:t_init] = fill(+Inf, nr, nk)
  out["GRID"][:t_init] = fill(+Inf, nr, nk)
  out["GRBS"][:t_init] = fill(+Inf, nr, nk)
  out["GRBU"][:t_init] = fill(+Inf, nr, nk)
  
  # t_run
  out["ETGS"][:t_run] = fill(+Inf, nr, nk)
  out["PLCP"][:t_run] = fill(+Inf, nr, nk)
  out["GRID"][:t_run] = fill(+Inf, nr, nk)
  out["GRBS"][:t_run] = fill(+Inf, nr, nk)
  out["GRBU"][:t_run] = fill(+Inf, nr, nk)

  # t_primal: primal recovery
  out["ETGS"][:t_primal] = fill(+Inf, nr, nk)
  out["PLCP"][:t_primal] = fill(+Inf, nr, nk)
  out["GRID"][:t_primal] = fill(+Inf, nr, nk)
  out["GRBS"][:t_primal] = fill(0.0, nr, nk)
  out["GRBU"][:t_primal] = fill(0.0, nr, nk)

  # t_sort
  out["ALL"][:t_sort] = +Inf # full sort time
  out["ALL"][:t_psort] = fill(+Inf, nk) # partial sort time
  out["GRBU"][:t_psort] = fill(0.0, nr, nk) # partial sort time
  
  # total infeas (sup-norm): tks and order
  out["ETGS"][:infeas] = fill(+Inf, nr, nk) # total sup infeasibility
  out["PLCP"][:infeas] = fill(+Inf, nr, nk) # total sup infeasibility
  out["GRID"][:infeas] = fill(+Inf, nr, nk) # total sup infeasibility
  out["GRBS"][:infeas] = fill(+Inf, nr, nk) # total sup infeasibility
  out["GRBU"][:infeas] = fill(+Inf, nr, nk) # total sup infeasibility
  
  # tks infeasibility
  out["ETGS"][:tks] = fill(+Inf, nr, nk) # method 1 violation of Mk(xbar) <= r
  out["PLCP"][:tks] = fill(+Inf, nr, nk) # method 2 violation of Mk(xbar) <= r
  out["GRID"][:tks] = fill(+Inf, nr, nk) # method 3 violation of Mk(xbar) <= r
  out["GRBS"][:tks] = fill(+Inf, nr, nk) # method 4 violation of Mk(xbar) <= r
  out["GRBU"][:tks] = fill(+Inf, nr, nk) # method 5 violation of Mk(xbar) <= r

  # order infeasibility
  out["ETGS"][:ord] = fill(+Inf, nr, nk) # method 1 violation of x_{i} >= x_{i+1}
  out["PLCP"][:ord] = fill(+Inf, nr, nk) # method 2 violation of x_{i} >= x_{i+1}
  out["GRID"][:ord] = fill(+Inf, nr, nk) # method 3 violation of x_{i} >= x_{i+1}
  out["GRBS"][:ord] = fill(+Inf, nr, nk) # method 4 violation of x_{i} >= x_{i+1}
  out["GRBU"][:ord] = fill(+Inf, nr, nk) # method 5 violation of x_{i} >= x_{[k]} for i in kappa and x_{i} <= x_{[k]} for i not in kappa
  
  # active constraints
  out["GRBS"][:nactivecon] = fill(0, nr, nk)
  out["GRBU"][:nactivecon] = fill(0, nr, nk)
  
  # obj
  out["ETGS"][:obj] = fill(+Inf, nr, nk) # obj value 1
  out["PLCP"][:obj] = fill(+Inf, nr, nk) # obj value 2
  out["GRID"][:obj] = fill(+Inf, nr, nk) # obj value 3
  out["GRBS"][:obj] = fill(+Inf, nr, nk) # obj value 4
  out["GRBU"][:obj] = fill(+Inf, nr, nk) # obj value 5
  
  # nit
  out["ETGS"][:nit] = zeros(Int64, nr, nk) # number of iterations 1
  out["PLCP"][:nit] = zeros(Int64, nr, nk) # number of iterations 2
  out["GRID"][:nit] = zeros(Int64, nr, nk) # number of iterations 3
  out["GRBS"][:nit] = zeros(Int64, nr, nk) # number of iterations 4
  out["GRBU"][:nit] = zeros(Int64, nr, nk) # number of iterations 5

  # case
  out["ETGS"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  out["PLCP"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  out["GRID"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  
  # bestfeas
  out["ETGS"][:bestfeas] = zeros(Bool, nr, nk) # does 1 satisfy Mk(xbar) <= r and have lowest objval
  out["PLCP"][:bestfeas] = zeros(Bool, nr, nk) # does 2 satisfy Mk(xbar) <= r and have lowest objval
  out["GRID"][:bestfeas] = zeros(Bool, nr, nk) # does 3 satisfy Mk(xbar) <= r and have lowest objval
  out["GRBS"][:bestfeas] = zeros(Bool, nr, nk) # does 4 satisfy Mk(xbar) <= r and have lowest objval
  out["GRBU"][:bestfeas] = zeros(Bool, nr, nk) # does 5 satisfy Mk(xbar) <= r and have lowest objval
  
  # other
  out["ALL"][:obj_12345] = Array{Vector{Int64}}(undef, nr, nk) # best objective value
  out["ALL"][:active] = zeros(Bool, nr, nk) # is Mk(x0) > r?
  out["ALL"][:x0sort] = zeros(DT, n) # x0sort
  
  # solution
  out["ETGS"][:xbarsort] = zeros(DT, n)
  out["PLCP"][:xbarsort] = zeros(DT, n)
  out["GRID"][:xbarsort] = zeros(DT, n)
  out["GRBS"][:xbarsort] = zeros(Float64, n)
  out["GRBU"][:xbarsort] = zeros(Float64, n)
  return out
end
function default_piv_options()
  D = Dict()
  D[:maxt] = 10_000
  D["ETGS"] = 1
  D["PLCP"] = 2
  D["GRID"] = 3
  D["GRBS"] = 4
  D["GRBU"] = 5
  D[1] = "ETGS"
  D[2] = "PLCP"
  D[3] = "GRID"
  D[4] = "GRBS"
  D[5] = "GRBU"
  return D
end

function writeout_piv(D::Dict, n::Integer, DT::DataType, repi::Integer, datapath::String)
  if isequal(D[:name], "ETGS")
    id = 1
  elseif isequal(D[:name], "PLCP")
    id = 2
  elseif isequal(D[:name], "GRID")
    id = 3
  else
    throw()
  end
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_k0_$repi.csv", D[:k0])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_k1_$repi.csv", D[:k1])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_init_$repi.csv", D[:t_init])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_run_$repi.csv", D[:t_run])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_primal_$repi.csv", D[:t_primal])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_total_$repi.csv", D[:t_init] + D[:t_run] + D[:t_primal])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_tks_$repi.csv", D[:tks])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_ord_$repi.csv", D[:ord])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_infeas_$repi.csv", D[:infeas])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_obj_$repi.csv", D[:obj])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_nit_$repi.csv", D[:nit])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_case_$repi.csv", D[:case])
end


function writeout_grb(D::Dict, n::Integer, DT::DataType, repi::Integer, datapath::String)
  if isequal(D[:name], "GRBS")
    id = 4
  elseif isequal(D[:name], "GRBU")
    id = 5
  else
    throw()
  end
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_k0_$repi.csv", D[:k0])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_k1_$repi.csv", D[:k1])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_init_$repi.csv", D[:t_init])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_run_$repi.csv", D[:t_run])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_primal_$repi.csv", D[:t_primal])
  # writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_total_$repi.csv", D[:t_init] + D[:t_run] + D[:t_primal])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_total_$repi.csv", D[:t_run] + D[:t_primal])
  if id == 5
    writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_psort_$repi.csv", D[:t_psort])
    # writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_total_$repi.csv", D[:t_init] + D[:t_run] + D[:t_primal] + D[:t_psort])
    writedlm(datapath*"/$(n)$(DT)/out_$(id)_t_total_$repi.csv", D[:t_run] + D[:t_primal])
  end
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_tks_$repi.csv", D[:tks])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_ord_$repi.csv", D[:ord])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_infeas_$repi.csv", D[:infeas])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_obj_$repi.csv", D[:obj])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_nit_$repi.csv", D[:nit])
  writedlm(datapath*"/$(n)$(DT)/out_$(id)_nactivecon_$repi.csv", D[:nactivecon])
end

"""
test_tks_timing(n::Integer, DT::DataType, nrep::Integer)
- n: problem dimension
- DT: data type
- nrep: number of replications
- rlevel: vector of tau_r
- klevel: vector of tau_k
- datapath: path for saving output
- maxn_gurobi: largest n for gurobi
- maxn_grid: largest n for grid
- expers: (r,k) tuples for subset of rlevel × klevel experiments
"""
function test_tks_timing(n::Integer, DT::DataType, nrep::Integer,
  rlevel::Vector, klevel::Vector,
  datapath::String,
  maxn_gurobi::Integer=100_000,
  maxn_grid::Integer=100_000,
  expers=nothing,
  rep_offset=0
)
  # setup
  nr = length(rlevel)
  nk = length(klevel)
  sig = DT(100) # x0_i ∈ [-sigma, +sigma]
  alg_names = ["ETGS", "PLCP", "GRID", "GRBS", "GRBU"]
  out = initialize_output(n, DT, nr, nk)
  
  # alg options
  piv_options = default_piv_options()
  piv_options[:maxt] = 10_000
  x0prepop = true # prepopulate xbarsort .= x0sort
  hist = false # don't record history
  verb = false # no console output
  
  # grb options
  grb_options = default_proj_options()
  grb_options[:maxtime] = 10_000
  @assert(maxn_gurobi <= maxn_grid) #! for `obj` vector ordering in determining `bestfeas`

  # test output
  println(datapath * "$n$DT/")
  mkpath(datapath * "/$n$DT/")
  writedlm(datapath*"/$(n)$(DT)/rlevel.csv", rlevel)
  writedlm(datapath*"/$(n)$(DT)/klevel.csv", klevel)
  writedlm(datapath*"/$(n)$(DT)/maxn_gurobi.csv", maxn_gurobi)
  writedlm(datapath*"/$(n)$(DT)/maxn_grid.csv", maxn_grid)
  
  sim_repi(repi::Int64, expers) = begin # modify `out` in parallel
    println("rep progress: $(round((repi-rep_offset) / nrep, digits=2))")
    total = length(klevel) * length(rlevel)
    remaining = total
    completed = 0
    printpct = 5
    flag = 0

    # grb environment
    begin
      env_p = Ref{Ptr{Cvoid}}()
      logfile = C_NULL
      error = GRBloadenv(env_p, logfile)
      @assert error == 0

      env = env_p[]

      error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0)
      @assert error == 0
      
      # parameters
      error = GRBsetdblparam(env, "TimeLimit", grb_options[:maxtime])
      @assert error == 0
      error = GRBsetintparam(env, "Method", Cint(2)) # barrier
      @assert error == 0
      error = GRBsetintparam(env, "BarIterLimit", grb_options[:maxiter])
      @assert error == 0
      error = GRBsetdblparam(env, "FeasibilityTol", grb_options[:pinfeas_tol])
      @assert error == 0
      error = GRBsetdblparam(env, "OptimalityTol", grb_options[:dinfeas_tol])
      @assert error == 0
      error = GRBsetdblparam(env, "BarConvTol", grb_options[:relgap_tol])
      @assert error == 0
      error = GRBsetintparam(env, "Threads", grb_options[:nthreads])
      @assert error == 0
    end
    
    # generate x0
    if DT <: Rational
      x0 = DT.(rand(0:sig, n)) .// sig
    elseif DT <: AbstractFloat
      x0 = rand(DT, n)
    end
    xbarsort = similar(x0);
    out["ALL"][:t_sort] = @elapsed begin
      x0sort = sort(x0, rev=true)
      out["ALL"][:x0sort] .= x0sort
    end
    xbarsort .= x0sort;
    writedlm(datapath*"/$(n)$(DT)/t_sort.csv", out["ALL"][:t_sort])

    # loop over k and r 
    for ki in eachindex(klevel)
      k = Int64(min(max(2, ceil(n*klevel[ki])),n-1))
      tks = sum(out["ALL"][:x0sort][1:k])
      out["ALL"][:t_psort][ki] = @elapsed begin
        partialsortsig = partialsortperm(x0, 1:k, rev=true)
      end
      writedlm(datapath*"/$(n)$(DT)/t_psort.csv", out["ALL"][:t_psort])
      
      for ri in eachindex(rlevel)
        r = DT(rlevel[ri]) * tks
        if !isnothing(expers)
          if (rlevel[ri], klevel[ki]) ∉ expers
            println("  --> skipping: $((rlevel[ri], klevel[ki]))")
            continue
          end
        end
        out["ALL"][:active][ri,ki] = sum(out["ALL"][:x0sort][1:k]) > r

        # ETGS
        xbarsort .= x0sort
        @views sp, k0k1, nit, t, case = project_topksum_esgs_experiment!(
           xbarsort, x0sort, r, k, out["ALL"][:active][ri,ki], x0prepop, verb, hist
        )
        D = out["ETGS"]
        D[:xbarsort] .= xbarsort
        D[:name] = "ETGS"
        D[:k0][ri,ki] = k0k1[1]
        D[:k1][ri,ki] = k0k1[2]
        D[:t_init][ri,ki] = t[1]
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3]
        D[:tks][ri,ki] = max(sum(D[:xbarsort][1:k]) - r, 0.0) #! is sorted
        D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
        D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
        D[:obj][ri,ki] = 0.5 * norm(D[:xbarsort] .- out["ALL"][:x0sort], 2)^2
        D[:nit][ri,ki] = nit
        D[:case][ri,ki] = case
        writeout_piv(D, n, DT, repi, datapath)

        # PLCP
        xbarsort .= x0sort
        @views sp, ab, nit, t, case = project_topksum_plcp_experiment!(
           xbarsort, x0sort, r, k, out["ALL"][:active][ri,ki], x0prepop, verb, hist
        )
        D = out["PLCP"]
        D[:xbarsort] .= xbarsort
        # k0k1 = get_k0k1(sort(D[:xbarsort], rev=true), k) #! soln not sorted but still can recover indices
        a = ab[1]
        b = ab[2]
        if a >=0 && b > 0
          a = max(a-1, 0)
          b = min(b+1, n)
        elseif case == 2
          a = k-1
          b = k
        end
        D[:name] = "PLCP"
        D[:k0][ri,ki] = a
        D[:k1][ri,ki] = b
        D[:t_init][ri,ki] = t[1]
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3]
        D[:tks][ri,ki] = max(sum(D[:xbarsort][1:k]) - r, 0.0) #! is sorted
        D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
        D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
        D[:obj][ri,ki] = 0.5 * norm(D[:xbarsort] .- out["ALL"][:x0sort], 2)^2
        D[:nit][ri,ki] = nit
        D[:case][ri,ki] = case
        writeout_piv(D, n, DT, repi, datapath)

        # GRID
        if n <= maxn_grid
          xbarsort .= x0sort
          @views sp, k0k1, nit, t, case = project_topksum_grid_experiment!(
           xbarsort, x0sort, r, k, out["ALL"][:active][ri,ki], x0prepop, verb, hist, piv_options[:maxt]
          )
          D = out["GRID"]
          D[:xbarsort] .= xbarsort
          D[:name] = "GRID"
          D[:k0][ri,ki] = k0k1[1]
          D[:k1][ri,ki] = k0k1[2]
          D[:t_init][ri,ki] = t[1]
          D[:t_run][ri,ki] = t[2]
          D[:t_primal][ri,ki] = t[3]
          D[:tks][ri,ki] = max(sum(sort(D[:xbarsort], rev=true)[1:k]) - r, 0.0) #! not necessarily sorted
          D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
          D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
          D[:obj][ri,ki] = 0.5 * norm(D[:xbarsort] .- out["ALL"][:x0sort], 2)^2
          D[:nit][ri,ki] = nit
          D[:case][ri,ki] = case
          writeout_piv(D, n, DT, repi, datapath)
        end
          
        # GUROBI
        if n <= maxn_gurobi
          # GRB sort
          res = project_topksum_grbs(out["ALL"][:x0sort], r, k, env, grb_options)
          D = out["GRBS"]
          D[:name] = "GRBS"
          D[:xbarsort] .= res[:x]
          k0k1 = get_k0k1(sort(D[:xbarsort], rev=true), k)
          D[:k0][ri,ki] = k0k1[1]
          D[:k1][ri,ki] = k0k1[2]
          D[:t_init][ri,ki] = res[:inittime]
          D[:t_run][ri,ki] = res[:walltime]
          D[:t_primal][ri,ki] = 0.0
          D[:tks][ri,ki] = max(sum(sort(D[:xbarsort], rev=true)[1:k]) - r, 0.0) #! not necessarily sorted
          D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
          D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
          D[:obj][ri,ki] = 0.5 * norm(D[:xbarsort] .- out["ALL"][:x0sort], 2)^2
          D[:nit][ri,ki] = res[:iter]
          D[:nactivecon][ri,ki] = res[:nactivecon]
          writeout_grb(D, n, DT, repi, datapath)
          
          # GRB unsort
          res = project_topksum_grbu(out["ALL"][:x0sort], r, k, env, grb_options)
          D = out["GRBU"]
          D[:name] = "GRBU"
          D[:xbarsort] .= res[:xsort]
          k0k1 = get_k0k1(sort(D[:xbarsort], rev=true), k)
          D[:k0][ri,ki] = k0k1[1]
          D[:k1][ri,ki] = k0k1[2]
          D[:t_init][ri,ki] = res[:inittime]
          D[:t_run][ri,ki] = res[:walltime]
          D[:t_primal][ri,ki] = 0.0
          D[:t_psort][ri,ki] = res[:psort_time]
          D[:tks][ri,ki] = max(sum(sort(D[:xbarsort], rev=true)[1:k]) - r, 0.0) #! not necessarily sorted
          kappa = res[:kappa]
          _k_ = res[:_k_]
          D[:ord][ri,ki] = max(
            max(0.0, maximum(res[:x][_k_] .- res[:x][kappa])),
            max(0.0, maximum(res[:x][setdiff(collect(1:n), kappa)] .- res[:x][_k_])),
          )
          D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
          D[:obj][ri,ki] = 0.5 * norm(D[:xbarsort] .- out["ALL"][:x0sort], 2)^2
          D[:nit][ri,ki] = res[:iter]
          D[:nactivecon][ri,ki] = res[:nactivecon]
          writeout_grb(D, n, DT, repi, datapath)
        end

        # method comparison
        feas = [
          out["ETGS"][:infeas][ri,ki]<=5eps(),
          out["PLCP"][:infeas][ri,ki]<=5eps(),
          out["GRID"][:infeas][ri,ki]<=5eps(),
          out["GRBS"][:infeas][ri,ki]<=5eps(),
          out["GRBU"][:infeas][ri,ki]<=5eps(),
        ]
        objs = [
          out["ETGS"][:obj][ri,ki],
          out["PLCP"][:obj][ri,ki],
          out["GRID"][:obj][ri,ki],
          out["GRBS"][:obj][ri,ki],
          out["GRBU"][:obj][ri,ki],
        ]
        out["ETGS"][:bestfeas][ri,ki] = feas[1] && prod(objs[1] .<= objs[feas].+5eps())
        out["PLCP"][:bestfeas][ri,ki] = feas[2] && prod(objs[2] .<= objs[feas].+5eps())
        out["GRID"][:bestfeas][ri,ki] = feas[3] && prod(objs[3] .<= objs[feas].+5eps())
        out["GRBS"][:bestfeas][ri,ki] = feas[4] && prod(objs[4] .<= objs[feas].+5eps())
        out["GRBU"][:bestfeas][ri,ki] = feas[5] && prod(objs[5] .<= objs[feas].+5eps())
        
        writedlm(datapath*"/$(n)$(DT)/out_1_bestfeas_$(repi).csv", out["ETGS"][:bestfeas])
        writedlm(datapath*"/$(n)$(DT)/out_2_bestfeas_$(repi).csv", out["PLCP"][:bestfeas])
        writedlm(datapath*"/$(n)$(DT)/out_3_bestfeas_$(repi).csv", out["GRID"][:bestfeas])
        writedlm(datapath*"/$(n)$(DT)/out_4_bestfeas_$(repi).csv", out["GRBS"][:bestfeas])
        writedlm(datapath*"/$(n)$(DT)/out_5_bestfeas_$(repi).csv", out["GRBU"][:bestfeas])
        writedlm(datapath*"/$(n)$(DT)/out_active_$(repi).csv", out["ALL"][:active])

        # check
        diff_12 = norm(out["ETGS"][:xbarsort] .- out["PLCP"][:xbarsort], Inf)
        diff_13 = norm(out["ETGS"][:xbarsort] .- out["GRID"][:xbarsort], Inf)
        diff_14 = norm(out["ETGS"][:xbarsort] .- out["GRBS"][:xbarsort], Inf)
        diff_15 = norm(out["ETGS"][:xbarsort] .- out["GRBU"][:xbarsort], Inf)
        diffs = [diff_12, diff_13, diff_14, diff_15]
        diffns = ["1v2", "1v3", "1v4", "1v5"]
        tol = 1e-3
        # sqrt(eps(DT)
        if argmax(diffs) == 1 && maximum(diffs) > tol
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in ETGS v PLCP)? tks infeas = $(
            round(max(0, sum(out["PLCP"][:xbarsort][1:k])-r),digits=5)
          ), diff_12 = $(
            round(diff_12,digits=8)
          )")
        elseif argmax(diffs) == 2 && maximum(diffs) > tol && n <= maxn_grid
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in ETGS v GRID)? tks infeas = $(
            round(max(0, sum(out["GRID"][:xbarsort][1:k])-r),digits=5)
          ), diff_13 = $(
            round(diff_13,digits=8)
          )")
        elseif argmax(diffs) == 3 && maximum(diffs) > tol && n <= maxn_gurobi
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in ETGS v GRBS)? tks infeas = $(
            round(max(0, sum(out["GRBS"][:xbarsort][1:k])-r),digits=5)
          ), diff_14 = $(
            round(diff_14,digits=8)
          )")
        elseif argmax(diffs) == 4 && maximum(diffs) > tol && n <= maxn_gurobi
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in ETGS v GRBU)? tks infeas = $(
            round(max(0, sum(out["GRBU"][:xbarsort][1:k])-r),digits=5)
          ), diff_15 = $(
            round(diff_15,digits=8)
          )")
        end
      end
    end
    GRBfreeenv(env_p[])
  end # sim_repi
  sr(i) = sim_repi(i + rep_offset, expers)
  pmap(sr, 1:nrep)
end

function load_tks_results(datapath::String, rep_lo::Integer, rep_hi::Integer,
  maxn_grid::Integer, maxn_gurobi::Integer, nlarge=false
)
  out = Dict()
  expers = readdir(datapath)
  expers = expers[.!occursin.(".csv", expers) .* .!occursin.(".tex", expers)]
  ns = [parse(Int64, match(r"(\d+)",x).match) for x in expers]
  exper_mask = ones(Bool, length(expers))
  if nlarge
    for i in eachindex(expers)
      exper = expers[i]
      files = readdir(datapath*exper*"/")
      rep_mask = [x[1] for x in split.(files, ".csv")]
      rep_mask = [x[end] for x in split.(rep_mask, "_")]
      rep_mask = [(isnothing(x) ? false : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
      if sum(rep_mask) == 0
        exper_mask[i] = false
      end
    end
  end

  for (exper,n) in zip(expers[exper_mask],ns[exper_mask])
    println("============================")
    println("exper = $exper")
    out[exper] = Dict()
    files = readdir(datapath*exper*"/")
    rep_mask = [x[1] for x in split.(files, ".csv")]
    rep_mask = [x[end] for x in split.(rep_mask, "_")]
    rep_mask = [(isnothing(x) ? true : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
    files = files[rep_mask]
    active = files[occursin.("_active", files)]
    bestfeas = files[occursin.("bestfeas", files)]
    case = files[occursin.("case", files)]
    infeas = files[occursin.("infeas", files)]
    k0 = files[occursin.("k0", files)]
    k1 = files[occursin.("k1", files)]
    tks = files[occursin.("tks", files)]
    nit = files[occursin.("_nit", files)]
    nactivecon = files[occursin.("nactivecon", files)]
    obj = files[occursin.("obj", files)]
    ord = files[occursin.("ord", files)]
    t_init = files[occursin.("t_init", files)]
    t_run = files[occursin.("t_run", files)]
    t_primal = files[occursin.("t_primal", files)]
    t_psort = files[occursin.("t_psort", files)]
    outnames = ["active", "bestfeas", "case", "infeas", "k0", "k1", "tks", "nit", "nactivecon", "obj", "ord", "t_init", "t_run", "t_primal", "t_psort"]
    outfiles = [active, bestfeas, case, infeas, k0, k1, tks, nit, nactivecon, obj, ord, t_init, t_run, t_primal, t_psort]
    klevel = files[occursin.("klevel", files)]
    rlevel = files[occursin.("rlevel", files)]
    # maxn_grid = files[occursin.("maxn_grid", files)]
    # maxn_gurobi = files[occursin.("maxn_gurobi", files)]
    nrep = maximum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_init)]) - minimum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_init)]) + 1
    nk = length(readdlm(datapath*exper*"/"*klevel[1]))
    nr = length(readdlm(datapath*exper*"/"*rlevel[1]))
    dims = (nrep,nr,nk)
    out[exper]["rlevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in readdlm(datapath*exper*"/"*rlevel[1])])
    out[exper]["klevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in readdlm(datapath*exper*"/"*klevel[1])])
    for algid in 1:5
      # println("----------------------------")
      # println("alg = $algid")
      if algid == 3 && maxn_grid < n
        continue
      end
      if (algid == 4 || algid == 5) && maxn_gurobi < n
        continue
      end
      out[exper][algid] = Dict()
      for (on,of) in zip(outnames,outfiles)
        # println(on)
        if on == "active"
          x = [readdlm(datapath*exper*"/"*f) for f in of[occursin.("out_", of)]]
          out[exper][on] = y = hvncat(dims, true, vcat(x...)...) # y[i,:,:]==x[i]
        elseif on == "case" && (algid == 4 || algid == 5)
          continue
        elseif (on == "nactivecon" || on == "t_psort") && (algid != 5)
          continue
        else
          x = [readdlm(datapath*exper*"/"*f) for f in of[occursin.("out_$(algid)_", of)]]
          out[exper][algid][on] = y = hvncat(dims, true, vcat(x...)...) # y[i,:,:]==x[i]
        end
        for i in 1:nrep
          @assert(y[i,:,:]==x[i])
        end
      end
      if algid == 4 || algid == 5
        out[exper][algid]["t_total"] = out[exper][algid]["t_run"]
      else
        out[exper][algid]["t_total"] = out[exper][algid]["t_init"] .+ out[exper][algid]["t_run"] .+ out[exper][algid]["t_primal"]
      end
    end
    println("-> $exper loaded")
  end
  out["sort_time"] = readdlm(datapath*"sort_time.csv")
  out["partial_sort_time"] = readdlm(datapath*"partial_sort_time_1pct.csv")
  # out["maxn_grid"] = parse(Int64,split(readdlm(datapath*"10Float64/maxn_grid.csv")[1],"//")[1])
  # out["maxn_gurobi"] = parse(Int64,split(readdlm(datapath*"10Float64/maxn_gurobi.csv")[1],"//")[1])
  out["rlevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in readdlm(datapath*"10Float64/rlevel.csv")])
  out["klevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in readdlm(datapath*"10Float64/klevel.csv")])
  return out
end

function calc_bestfeasunsorted!(out::Dict)
  for key in keys(out)
    println("exper = $key")
    if !isnothing(match(r"(\d+)", key))
      for alg in 1:5
        if haskey(out[key], alg)
          nrep = size(out[key][alg]["tks"], 1)
          bestobj = copy(out[key][alg]["obj"])
          for i in 1:5
            if haskey(out[key], i)
              bestobj .= min.(bestobj, out[key][i]["obj"] .+ 1e10 .* (out[key][i]["tks"] .> 5eps()))
            end
          end
          out[key][alg]["bestfeasunsorted"] = (
            (out[key][alg]["tks"] .<= 5eps()) .* (out[key][alg]["obj"] .<= bestobj .+ 5eps())
          )
        end
      end
    end
  end
end
