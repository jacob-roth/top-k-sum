function default_proj_options()
  options = Dict()
  options[:maxtime] = 10_000 # seconds
  options[:method] = 2 # barrier
  options[:maxiter] = 10_000
  options[:pinfeas_tol] = 1e-9
  options[:dinfeas_tol] = 1e-9
  options[:relgap_tol] = 1e-9
  options[:nthreads] = 8
  return options
end

function project_topksum_grbs(
  x0sort::Vector, r::Real, k::Integer,
  env::Ptr{Nothing},
  options::Dict=default_proj_options()
)
  start_time = time()
  
  #
  # initialize
  #

  # dims
  n = length(x0sort)
  n_Cint = Cint(n)
  n_Csize_t = Base.Csize_t(n)
  xidx = 1:n

  #
  # A matrix
  #

  indk = hcat(ones(k)', spzeros(1, n-k))
  B = spdiagm(0=>ones(n), 1=>-ones(n-1))[1:n-1,:]
  A = [
  # x
    +indk./k # c1: indk' * x <= r
    B        # c3: Bx >= 0
  ]
  c1idx = 1:1
  c2idx = 2:n
  Acsr = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(A))));
  
  #
  # environment
  #

  # parameters
  error = GRBsetdblparam(env, "TimeLimit", options[:maxtime])
  @assert error == 0
  error = GRBsetintparam(env, "Method", Cint(2)) # barrier
  @assert error == 0
  error = GRBsetintparam(env, "BarIterLimit", options[:maxiter])
  @assert error == 0
  error = GRBsetdblparam(env, "FeasibilityTol", options[:pinfeas_tol])
  @assert error == 0
  error = GRBsetdblparam(env, "OptimalityTol", options[:dinfeas_tol])
  @assert error == 0
  error = GRBsetdblparam(env, "BarConvTol", options[:relgap_tol])
  @assert error == 0
  error = GRBsetintparam(env, "Threads", options[:nthreads])
  @assert error == 0
  
  #
  # model
  #

  model_p = Ref{Ptr{Cvoid}}()
  modelname = "proj"
  error = GRBnewmodel(env, model_p, modelname, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
  @assert error == 0
  model = model_p[]
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)
  @assert error == 0

  #
  # bounds
  #

  x_lb = fill(typemin(Cdouble), n)
  x_ub = fill(typemax(Cdouble), n)
  
  #
  # linear objective
  #

  x_lobj = zeros(Cdouble, n)
  x_lobj .= -x0sort
  
  #
  # variables
  #

  error = GRBXaddvars(
      model,      # : model
      n_Csize_t,  # : numvars
      0,          # : numnz
      C_NULL,     # : *vbeg
      C_NULL,     # : *vind
      C_NULL,     # : *vval
      x_lobj,     # : *obj
      x_lb,       # : *lb
      x_ub,       # : *ub
      C_NULL,     # : *vtype
      C_NULL      # : **varnames
  )
  @assert error == 0

  #
  # quadratic objective
  #

  qrow, qcol, qval = findnz(sparse(spdiagm(ones(n))))
  error = GRBaddqpterms(
    model,             # GRBmodel : *model
    Cint(n),           # int      : numqnz
    Cint.(qrow.-1),    # int      : *qrow
    Cint.(qcol.-1),    # int      : *qcol
    0.5Cdouble.(qval), # double   : *qval #! why 0.5???
  )
  @assert error == 0

  #
  # constraints
  #

  consense = fill(Cchar(0), Acsr.m)
  consense[c1idx] .= Cchar(GRB_LESS_EQUAL)
  consense[c2idx] .= Cchar(GRB_GREATER_EQUAL)
  conname = fill(C_NULL, Acsr.m)
  rhs = [
    +Cdouble(r/k); # c1
    zeros(Cdouble, n-1); # c2
  ]
  error = GRBXaddconstrs(
    model,                         # GRBmodel   : *model
    Cint(Acsr.m),                  # int        : numconstrs
    Base.Csize_t(nnz(Acsr)),       # size_t     : numnz
    Base.Csize_t.(Acsr.rowptr.-1), # size_t     : *cbeg
    Cint.(Acsr.colval.-1),         # int        : *cind
    Cdouble.(Acsr.nzval),          # double     : *cval
    consense,                      # char       : sense
    rhs,                           # double     : *rhs
    conname,                       # const char : **constrname
  )
  @assert error == 0
  
  #
  # update solver parameters
  #

  # update
  GRBupdatemodel(model)
  time_init = time()-start_time
  # println("model initialization = $(time()-start_time)")

  #
  # solve
  #

  solve_t = @elapsed error = GRBoptimize(model)
  # println("model solve = $(solve_t)")
  
  #
  # recover solution
  #

  res = Dict();
  NumVars = Ref{Cint}();
  NumConstrs = Ref{Cint}();
  IterCount = Ref{Cint}(); # simplex iters
  BarIterCount = Ref{Cint}(); # barrier iters
  ObjVal = Ref{Cdouble}(); # primal obj
  ObjBound = Ref{Cdouble}(); # dual obj
  RunTime = Ref{Cdouble}(); # solver wall-clock time
  x = zeros(Cdouble, n); # x
  Status = Ref{Cint}(); # solved status
  Slack = zeros(Cdouble, Acsr.m);
  
  error = GRBgetintattr(model, "NumVars", NumVars);
  @assert error == 0
  error = GRBgetintattr(model, "NumConstrs", NumConstrs);
  @assert error == 0
  error = GRBgetintattr(model, "Status", Status);
  @assert error == 0
  error = GRBgetintattr(model, "BarIterCount", BarIterCount);
  @assert error == 0
  error = GRBgetdblattr(model, "ObjVal", ObjVal);
  @assert error == 0
  error = GRBgetdblattr(model, "RunTime", RunTime);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, x);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_SLACK, 0, Acsr.m, Slack);
  @assert error == 0
  
  # populate result
  res[:numvar] = NumVars[]
  res[:numcon] = NumConstrs[]
  res[:walltime] = RunTime[]
  res[:inittime] = time_init
  method = 2
  if method == 2
    res[:iter] = BarIterCount[]
  else
    res[:iter] = IterCount[]
  end
  res[:pobj] = ObjVal[]
  res[:x] = x[xidx]
  res[:status] = (Status[] == 2 ? 1 : 0)
  res[:nactivecon] = sum(abs.(Slack) .<= 10options[:pinfeas_tol])
  
  #
  # free model
  #
  
  @assert error == 0
  error = GRBfreemodel(model)
  @assert error == 0
  return res
end


function project_topksum_grbu(
  x0::Vector, r::Real, k::Integer,
  env::Ptr{Nothing},
  options::Dict=default_proj_options()
)
  # not including partial sort time; show it explicitly as a line in the table
  psort_time = @elapsed begin 
    sig = sortperm(x0, alg=PartialQuickSort(1:k), rev=true);
    kappa = sig[1:k]
  end
  _k_ = kappa[end]

  # start timer
  start_time = time()
  
  #
  # initialize
  #

  # dims
  n = length(x0)
  n_Cint = Cint(n)
  n_Csize_t = Base.Csize_t(n)
  xidx = 1:n

  #
  # A matrix
  #

  indkappa = spzeros(1, n)
  indkappa[kappa] .= 1
  col = [1:_k_-1; _k_.*ones(Int64, _k_-1);  _k_.*ones(Int64, n-_k_); _k_+1:n]
  nelptr = [1; ones(Int64, _k_-1); n-1; ones(Int64, n-_k_)]
  colptr = cumsum(nelptr)
  rowval = [1:_k_-1; 1:_k_-1; _k_+1:n; _k_+1:n]
  nzval = [-ones(_k_-1); ones(_k_-1); -ones(n-_k_); ones(n-_k_)]
  A = SparseMatrixCSC(n, n, colptr, rowval, nzval)
  A = A[[1:_k_-1; _k_+1:n],:]
  B = [
    +indkappa./k # c1: (indkappa' * x)/k <= r/k
    A            # c2: Ax <= 0
  ]
  c1idx = 1:1
  c2idx = 2:n
  Bcsr = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(B))));

  #
  # environment
  #
  
  # parameters
  error = GRBsetdblparam(env, "TimeLimit", options[:maxtime])
  @assert error == 0
  error = GRBsetintparam(env, "Method", Cint(2)) # barrier
  @assert error == 0
  error = GRBsetintparam(env, "BarIterLimit", options[:maxiter])
  @assert error == 0
  error = GRBsetdblparam(env, "FeasibilityTol", options[:pinfeas_tol])
  @assert error == 0
  error = GRBsetdblparam(env, "OptimalityTol", options[:dinfeas_tol])
  @assert error == 0
  error = GRBsetdblparam(env, "BarConvTol", options[:relgap_tol])
  @assert error == 0
  error = GRBsetintparam(env, "Threads", options[:nthreads])
  @assert error == 0

  #
  # model
  #

  model_p = Ref{Ptr{Cvoid}}()
  modelname = "proj"
  error = GRBnewmodel(env, model_p, modelname, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
  @assert error == 0
  model = model_p[]
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)
  @assert error == 0

  #
  # bounds
  #

  x_lb = fill(typemin(Cdouble), n)
  x_ub = fill(typemax(Cdouble), n)
  
  #
  # linear objective
  #

  x_lobj = zeros(Cdouble, n)
  x_lobj .= -x0
  
  #
  # variables
  #

  error = GRBXaddvars(
      model,      # : model
      n_Csize_t,  # : numvars
      0,          # : numnz
      C_NULL,     # : *vbeg
      C_NULL,     # : *vind
      C_NULL,     # : *vval
      x_lobj,     # : *obj
      x_lb,       # : *lb
      x_ub,       # : *ub
      C_NULL,     # : *vtype
      C_NULL      # : **varnames
  )
  @assert error == 0

  #
  # quadratic objective
  #

  qrow, qcol, qval = findnz(sparse(spdiagm(ones(n))))
  error = GRBaddqpterms(
    model,             # GRBmodel : *model
    Cint(n),           # int      : numqnz
    Cint.(qrow.-1),    # int      : *qrow
    Cint.(qcol.-1),    # int      : *qcol
    0.5Cdouble.(qval), # double   : *qval #! why 0.5???
  )
  @assert error == 0

  #
  # constraints
  #

  consense = fill(Cchar(0), Bcsr.m)
  consense[c1idx] .= Cchar(GRB_LESS_EQUAL)
  consense[c2idx] .= Cchar(GRB_LESS_EQUAL)
  conname = fill(C_NULL, Bcsr.m)
  rhs = [
    +Cdouble(r/k); # c1
    zeros(Cdouble, n-1); # c3
  ]
  error = GRBXaddconstrs(
    model,                         # GRBmodel   : *model
    Cint(Bcsr.m),                  # int        : numconstrs
    Base.Csize_t(nnz(Bcsr)),       # size_t     : numnz
    Base.Csize_t.(Bcsr.rowptr.-1), # size_t     : *cbeg
    Cint.(Bcsr.colval.-1),         # int        : *cind
    Cdouble.(Bcsr.nzval),          # double     : *cval
    consense,                      # char       : sense
    rhs,                           # double     : *rhs
    conname,                       # const char : **constrname
  )
  @assert error == 0
  
  #
  # update solver parameters
  #

  # update
  GRBupdatemodel(model)
  time_init = time()-start_time
  # println("model initialization = $(time()-start_time)")

  #
  # solve
  #

  solve_t = @elapsed error = GRBoptimize(model)
  # println("model solve = $(solve_t)")
  
  #
  # recover solution
  #

  res = Dict();
  NumVars = Ref{Cint}();
  NumConstrs = Ref{Cint}();
  IterCount = Ref{Cint}(); # simplex iters
  BarIterCount = Ref{Cint}(); # barrier iters
  ObjVal = Ref{Cdouble}(); # primal obj
  ObjBound = Ref{Cdouble}(); # dual obj
  RunTime = Ref{Cdouble}(); # solver wall-clock time
  x = zeros(Cdouble, n); # x
  Slack = zeros(Cdouble, Bcsr.m); # x
  Status = Ref{Cint}(); # solved status
  
  error = GRBgetintattr(model, "NumVars", NumVars);
  @assert error == 0
  error = GRBgetintattr(model, "NumConstrs", NumConstrs);
  @assert error == 0
  error = GRBgetintattr(model, "Status", Status);
  @assert error == 0
  error = GRBgetintattr(model, "BarIterCount", BarIterCount);
  @assert error == 0
  error = GRBgetdblattr(model, "ObjVal", ObjVal);
  @assert error == 0
  error = GRBgetdblattr(model, "RunTime", RunTime);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, x);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_SLACK, 0, Bcsr.m, Slack);
  @assert error == 0
  
  # populate result
  res[:numvar] = NumVars[]
  res[:numcon] = NumConstrs[]
  res[:walltime] = RunTime[]
  res[:inittime] = time_init
  res[:psort_time] = psort_time
  method = 2
  if method == 2
    res[:iter] = BarIterCount[]
  else
    res[:iter] = IterCount[]
  end
  res[:pobj] = ObjVal[]
  res[:x] = x[xidx]
  res[:kappa] = kappa
  res[:_k_] = _k_
  res[:nactivecon] = sum(abs.(Slack) .<= 10options[:pinfeas_tol])
  res[:xsort] = sort(res[:x], rev=true)
  res[:status] = (Status[] == 2 ? 1 : 0)
  
  #
  # free model
  #
  
  @assert error == 0
  error = GRBfreemodel(model)
  @assert error == 0
  return res
end
