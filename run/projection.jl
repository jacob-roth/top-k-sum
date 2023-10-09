#=================================================
PLCP: detailed timing of Algorithm 1
  - xbarsort: preallocated output vector
  - x0sort: sorted input vector
  - r: scalar budget parameter
  - k: integer
  - active: binary for whether Tk(x0sort) > r
  - x0prepop: binary for whether xbarsort .= x0sort
  - verb: binary for verbosity level
  - hist: binary for history
=================================================#
function project_topksum_plcp_experiment!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti, active::Bool=true,
  x0prepop::Bool=false, verb::Bool=false, hist::Bool=false
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
time_init_1 = @elapsed begin
  # inactive
  n = length(x0sort)
  if !active
    case = 0
    if !x0prepop
      @simd for i in 1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
    if hist
      return 0, (-1, -1), 0, (0.0, 0.0, 0.0), [(0,0)], case, 0.0
    else
      return 0, (-1, -1), 0, (0.0, 0.0, 0.0), case, 0.0
    end
  end

  # initialize
  solved::Bool = false
  n = length(x0sort)
  t::Ti = 0
  s0::Tfr = sum(view(x0sort, 1:k))
  qk::Tfr = (k == n ? typemax(Tfr) : x0sort[k] - x0sort[k+1])
  lam::Tfr = 0
  case::Ti = 0
end # time_init_1

  if k == n
    if verb
      println("case k=n: solve directly")
    end
time_primal = @elapsed begin
    # xbarsort = x0sort - (s0 - r)/n
    lam = (s0 - r) / k
    case = -2
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lam
    end
end # time_primal
    if hist
      return 0, (0,n), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, (0,n), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  elseif k == 1
    if verb
      println("case k=1: solve directly")
    end
time_primal = @elapsed begin
    case = -1
    # xbarsort = min(x0sort, r)
    if x0prepop
      for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
        if x0sort[i] <= r
          break
        end
      end
    else
      @simd for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
      end
    end
end # time_primal
    if hist
      return 0, (0,findfirst(x0sort .< r)), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, (0,findfirst(x0sort .< r)), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  else
    # t = 0
    # evaluate mks @ x(lam_1): lam_1=q_k & z=0
time_primal = @elapsed begin
    case = 2
    lam = qk
    m = s0 - k * lam
    if m <= r
      lam = (s0 - r) / k
      if verb
        println("lam = $lam")
        println("m = s0 - k * lam = $s0 - $(k * lam) = $(s0 - k * lam)")
        println("r - m >= 0? r-m = $(r-m)")
        println("solved t=0!")
        println("lamstar=$lam")
      end
      solved = true
      if x0prepop
        @simd for i in 1:k
          @inbounds xbarsort[i] -= lam
        end
      else
        @simd for i in 1:k
          @inbounds xbarsort[i] = x0sort[i] - lam
        end
        @simd for i in k+1:n
          @inbounds xbarsort[i] = x0sort[i]
        end
      end
    end
end # time_primal
    if solved
      if hist
        return 0, (0, 0), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
      else
        return 0, (0, 0), 0, (time_init_1, 0.0, time_primal), case, lam
      end
    end
  end

time_init_2 = @elapsed begin
  # t = 1
  case = 1
  t = 1
  a = k
  b = k
  a_alpha = one(Ti)
  k_alpha = one(Ti)
  b_alpha = one(Ti)
  sigma::Tfr = 3 // 2
  lam = qk
  lam_a::Tfr = qk
  lam_b::Tfr = qk
  za::Tfr = -qk / 2
  zk::Tfr = -qk / 2 # solve 2zk - 1zk+1 + qk = 0
  zb::Tfr = -qk / 2
  minv_ak = (t + one(Tfr) - k_alpha) / (t + one(Tfr)) # == minv_ka
  minv_bk = one(Tfr) - minv_ak # == minv_kb
  minv_ab = (t + one(Tfr) - b_alpha) / (t + one(Tfr)) # == minv_ba
  minv_kk = ((t + one(Tfr) - k_alpha) * k_alpha ) / (t + one(Tfr))
  pInf = typemax(Tfr)
  nInf = typemin(Tfr)

  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (a, b)
  end
end # time_init_2

  # iterate
time_pivot = @elapsed begin
  while true

    if verb
      println("\n-------------------------")
      println("t = $t")
      println("(a,b) = ($a,$b)")
      println("(za, zk, zb)  = $((za, zk, zb))")
    end

    # get breakpoint
    lam_a = (a > 1 ? ((x0sort[a-1] - x0sort[a]) - za) / minv_ak : pInf) # (q[a-1] - za) / minv_ka
    if b+1 < n
      lam_b = ((x0sort[b+1] - x0sort[b+2]) - zb) / minv_bk
    else
      lam_b = pInf
    end
    lam = min(lam_a, lam_b)
    if verb
      println("lam_a = $lam_a")
      println("lam_b = $lam_b")
      if lam_a <= lam_b
        println("s = a")
      else
        println("s = b")
      end
      println("lam = $lam\n")
    end
    
    # check max-k-sum
    if lam != pInf
      m = s0 - k * lam + zk + lam * minv_kk
      if verb
        println("m = s0 - k * lam + zk + lam * minv_kk = $m")
        println("r-m = $(r-m)")
        println("zk = $zk, lam=$lam, minv_kk=$minv_kk")
        println("compare to m' = s0 - k*lam + (-qk + lam)/2 = $(s0 - k * lam + (-qk + lam)/2)")
        println("r-m' = $(r - ((s0 - k * lam + (-qk + lam)/2)))")
      end
    else
      m = nInf
    end
    if m <= r# + tol
      lam = (s0 - r + zk) / (k - minv_kk)
      if verb
        println("m <= r, so rootfinding")
        println("m=$m, r=$r")
        println("(s0 - r + zk) = $s0 - $r + $zk = $(s0 - r + zk)")
        println("k - minv_kk = $k - $minv_kk = $(k-minv_kk)")
        println("lamstar = $lam")
      end
      solved = true
      break
    end
    
    # update solution
    if lam_a <= lam_b
      # s = a-1
      if verb; println("za = ($za - ($(x0sort[a-1]) - $(x0sort[a]))) / $sigma"); end
      za = (za - (x0sort[a-1] - x0sort[a])) / sigma # q[a-1] = x0sort[a-1] - x0sort[a]
      zk = zk + za * minv_ak
      zb = zb + za * minv_ab
      if verb; println("zk = $zk + $za * $minv_ak"); end
      if verb; println("zb = $zb + $za * $minv_ab"); end
      a -= 1
      k_alpha += 1
      b_alpha += 1
    else
      # s = b+1
      if verb; println("zb = ($zb - ($(x0sort[b+1]) - $(x0sort[b+2]))) / $sigma"); end
      zb = (zb - (x0sort[b+1] - x0sort[b+2])) / sigma # q[b+1] = x0sort[b+1] - x0sort[b+2]
      za = za + zb * minv_ab
      zk = zk + zb * minv_bk
      if verb; println("za = $za + $zb * $minv_ab"); end
      if verb; println("za = $zk + $zb * $minv_bk"); end
      b += 1
      b_alpha += 1
    end

    # increment basis length
    t += 1 #! number of iterations == t-1

    if hist
      history[t] = (a, b)
    end

    # update M(alpha)⁻¹ quantities
    minv_ak = (t + one(Tfr) - k_alpha) / (t + one(Tfr))
    minv_bk = one(Tfr) - minv_ak
    minv_ab = (t + one(Tfr) - b_alpha) / (t + one(Tfr))
    minv_kk = ((t + one(Tfr) - k_alpha) * k_alpha ) / (t + one(Tfr))
    if Tfr<:Rational
      sigma = (t+2) // (t+1)
    else
      sigma = (t+2) / (t+1)
    end
    if verb
      println("t / (t+1) = $(t / (t + 1))")
      println("sigma=$sigma")
    end
  end
end # time_pivot

time_primal = @elapsed begin
  if solved
    if !x0prepop
      xbarsort .= x0sort
    end
    for i in 1:k
      xbarsort[i] -= lam
    end
    
    ypBtz!(xbarsort, x0sort, a, b, k_alpha, lam)
  end
end # time_primal

  # return
  if hist
    return 1, (a,b), t, (time_init_1+time_init_2, time_pivot, time_primal), hist, case, lam #! not t-1 bc haven't incremented t yet
  else
    return 1, (a,b), t, (time_init_1+time_init_2, time_pivot, time_primal), case, lam #! not t-1 bc shaven't incremented t yet
  end
  nothing
end

#=================================================
ESGS: detailed timing of Algorithm 3
  - xbarsort: preallocated output vector
  - x0sort: sorted input vector
  - r: scalar budget parameter
  - k: integer
  - active: binary for whether Tk(x0sort) > r
  - x0prepop: binary for whether xbarsort .= x0sort
  - verb: binary for verbosity level
  - hist: binary for history
=================================================#
function project_topksum_esgs_experiment!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, verb::Bool=false, hist::Bool=false,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
time_init_1 = @elapsed begin
  # inactive
  n = length(x0sort)
  if !active
    case = 0
    if !x0prepop
      @simd for i in 1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (0.0, 0.0, 0.0), [(0,0)], case, 0.0
    else
      return 0, get_k0k1(xbarsort, k), 0, (0.0, 0.0, 0.0), case, 0.0
    end
  end

  # initialize
  solved::Bool = false
  n = length(x0sort)
  s0::Tfr = sum(view(x0sort, 1:k))
  lam::Tfr = 0
  case::Ti = 0
end # time_init_1

  if k == n
    if verb
      println("case k=n: solve directly")
    end
time_primal = @elapsed begin
    case = -2
    # xbarsort = x0sort - (s0 - r)/n
    lam = (s0 - r) / k
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lam
    end
end # time_primal
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  elseif k == 1
    if verb
      println("case k=1: solve directly")
    end
time_primal = @elapsed begin
    case = -1
    # xbarsort = min(x0sort, r)
    if x0prepop
      for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
        if x0sort[i] <= r
          break
        end
      end
    else
      @simd for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
      end
    end
end # time_primal
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  end

time_init = @elapsed begin
  # preprocessing
  tol = zero(Tfr)
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    # tol = eps(Tfr)*maximum(x0sort)
    tol = eps(Tfr)*x0sort[1]
  end
  case = 1
  n::Ti = length(x0sort)
  k0::Ti = k-1
  k1::Ti = k
  kk0::Ti = 0
  k1k0::Ti = 0
  the = zero(Tfr) # theta
  lam = zero(Tfr) # lambda
  theplam = zero(Tfr) # theta + lambda
  kkt2 = false
  kkt5 = false
  if verb
    kkt1 = false
    kkt3 = false
    kkt4 = false
  end
  nit::Ti = 0
  maxnit::Ti = length(x0sort)
  sum1k0::Tfr = s0 - x0sort[k]
  sumk0p1k1::Tfr = x0sort[k]
  sum1k0r = zero(Tfr)
  rho = zero(Tfr)
  solved = false
  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (k0, k1)
  end
end # time_init

  # iterate
time_pivot = @elapsed begin
  while true

    # counter
    nit += 1

    # rho, theta, and lambda
    kk0 = k - k0
    k1k0 = k1 - k0
    rho = k0 * k1k0 + kk0^2
    sum1k0r = sum1k0 - r
    the = (k0 * sumk0p1k1 - kk0 * sum1k0r) / rho #! OK
    if k0 > 0
      theplam = (k * the + sum1k0r) / k0
    else
      theplam = (k * sumk0p1k1 + (k1-k) * sum1k0r) / rho
    end

    # kkt conditions
    kkt2 = (k0 == 0 ? true : x0sort[k0] > theplam-tol) #! OK
    kkt5 = (k1 == n ? true : the > x0sort[k1+1]-tol) #! OK
    
    # debug
    if verb
      kkt1 = lam > 0
      kkt3 = the + lam >= x0sort[k0+1]
      kkt4 = x0sort[k1] >= the
      println("----------------")
      println("k0=$k0, k1=$k1")
      println("the=$the, lam=$(theplam-the), rho=$rho")
      println("sum1k0r_rho=$(sum1k0r/rho), sumk0p1k1_rho=$(sumk0p1k1/rho)")
      println("kkt1=$kkt1, val=$(theplam-the)")
      if k0 > 0
        println("kkt2=$kkt2, val=$(k0==0 ? Inf : x0sort[k0]-(theplam))", ", $(x0sort[k0]) > $(theplam))")
      end
      println("kkt3=$kkt3, val=$((theplam) - x0sort[k0+1])")
      println("kkt4=$kkt4, val=$(x0sort[k1] - the)")
      if k1 < n
        println("kkt5=$kkt5, val=$(the - (k1==n ? -Inf : x0sort[k1+1]))")
      end
    end

    # move k0
    if (kkt2 && kkt5) || (k0 == 0 && k1 == n)
      solved = true
      break
    elseif kkt2 && (k1 < n)
      k1 += 1
      sumk0p1k1 += x0sort[k1]
      if verb
        printstyled("k1 <- k1+1\n"; color=:green)
      end
    elseif !kkt2 && (k0 > 0)
      sum1k0 -= x0sort[k0]
      sumk0p1k1 += (k0 > 0 ? x0sort[k0] : 0)
      k0 -= 1
      if verb
        printstyled("k0 <- k0-1\n"; color=:blue)
      end
    end

    # history
    if hist
      history[nit+1] = (k0, k1)
    end
  end
end # time_pivot

  # stop
time_primal = @elapsed begin
  if solved
    lam = theplam - the
    if x0prepop
      @simd for i in 1:k0
        @inbounds xbarsort[i] -= lam
      end
    else
      @simd for i in 1:k0
        @inbounds xbarsort[i] = x0sort[i] - lam
      end
    end
    @simd for i in k0+1:k1
      @inbounds xbarsort[i] = the
    end
    if !x0prepop
      @simd for i in k1+1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
  end
end # time_primal

  if solved
    if !hist
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), case, lam
    else
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), history, case, lam
    end
  elseif nit >= maxnit
    if !hist
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), case, lam
    else
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), history, case, lam
    end
  end
end

#=================================================
GRID: detailed timing of GRID algorithm
  - xbarsort: preallocated output vector
  - x0sort: sorted input vector
  - r: scalar budget parameter
  - k: integer
  - active: binary for whether Tk(x0sort) > r
  - x0prepop: binary for whether xbarsort .= x0sort
  - verb: binary for verbosity level
  - hist: binary for history
=================================================#
function project_topksum_grid_experiment!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, verb::Bool=false, hist::Bool=false, maxt::Real=10_000,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
t_start = time()
time_init_1 = @elapsed begin
  # inactive
  n = length(x0sort)
  if !active
    case = 0
    if !x0prepop
      @simd for i in 1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
    if hist
      return 0, (-1, -1), 0, (0.0, 0.0, 0.0), [(0,0)], case, 0.0
    else
      return 0, (-1, -1), 0, (0.0, 0.0, 0.0), case, 0.0
    end
  end

  # initialize
  solved::Bool = false
  n = length(x0sort)
  s0::Tfr = sum(x0sort[1:k])
  lambda::Tfr = 0
  case::Ti = 0
end # time_init_1

  if k == n
    if verb
      println("case k=n: solve directly")
    end
time_primal = @elapsed begin
    case = -2
    # xbarsort = x0sort - (s0 - r)/n
    lambda = (s0 - r) / k
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lambda
    end
end # time_primal
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lambda
    else
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), case, lambda
    end
  elseif k == 1
    if verb
      println("case k=1: solve directly")
    end
time_primal = @elapsed begin
    case = -1
    # xbarsort = min(x0sort, r)
    if x0prepop
      for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
        if x0sort[i] <= r
          break
        end
      end
    else
      @simd for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
      end
    end
end # time_primal
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, get_k0k1(xbarsort, k), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  end

time_init_1 = @elapsed begin
  #=
  preprocessing
  =#
  case = 0
  tol::Tfr = 0
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    tol = eps(Tfr)*x0sort[1]
  end
  n = length(x0sort)
  nit::Ti = 0
  
  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (k0, k1)
  end
end # time_init_1

  #= 
  case 1
  =#
  case = 1
  time_init_2 = @elapsed begin
  # initialize
  k0::Ti = k - one(Ti);
  flag::Ti = zero(Ti);
  k1::Ti = k;
  sum1k = s0;
  sum1k0::Tfr = sum1k - x0sort[k];
  sumk0p1k::Tfr = x0sort[k];
  sumk0p1k1::Tfr = zero(Tfr);
  sumk0p1k1 = sumk0p1k;
  rho::Ti = k0 * (k1 - k0) + (k - k0)^2;
  kk0::Ti = k - k0;
  k1k0::Ti = k1 - k0;
  sum1k0r_rho::Tfr = zero(Tfr); # (sum1k0 - r) / rho
  sumk0p1k1_rho::Tfr = zero(Tfr); # sumk0p1k1 / rho
  theta::Tfr = zero(Tfr);
  R::Tfr = zero(Tfr); # R := 
end # time_init_2
  
  # iterate
time_pivot = @elapsed begin
  while true
    if verb
      println("------------------------------")
      println("k0=$k0, k1=$k1")
    end
    if time() - t_start > maxt
      println("!TIMELIMIT!")
      case = -100
      if hist
        return 1, (k0, k1), nit, (time_init_1+time_init_2, time()-t_start, 0.0), hist, case, lambda
      else
        return 1, (k0, k1), nit, (time_init_1+time_init_2, time()-t_start, 0.0), case, lambda
      end
      break
    end
    nit += 1;
    kk0 = k - k0;
    k1k0 = k1 - k0;
    rho = k0 * k1k0 + kk0^2;
    
    sum1k0r_rho = (sum1k0 - r) / rho;
    sumk0p1k1_rho = sumk0p1k1 / rho;
    theta = k0 * sumk0p1k1_rho - kk0 * sum1k0r_rho;
    lambda = kk0 * sumk0p1k1_rho + k1k0 * sum1k0r_rho;
    
    if verb
      println("kkt1: ", (lambda > -tol)," ", lambda)
      if k0 > 0
        println("kkt2: ", (k0 == 0 ? true : x0sort[k0] > theta + lambda - tol)," $(x0sort[k0]) > $(theta + lambda - tol)")
      end
      println("kkt3: ", (theta + lambda >= x0sort[k0+1] - tol)," $(theta+lambda) > $(x0sort[k0+1])")
      println("kkt4: ", (x0sort[k1] >= theta - tol)," $(x0sort[k1]) >= $theta")
      if k1 < n
        println("kkt5: ", (k1+1 > n ? true : theta > x0sort[k1+1] - tol)," $theta > $(x0sort[k1+1])")
      end
    end

    if (k0 == 0) && (k1 == n)
      flag = 1
    elseif ((k0==k-1 && k1==k) && 
        (k0 == 0 ? true : x0sort[k0] > theta + lambda - tol) && #2
        (k1 == n ? true : theta > x0sort[k1+1] - tol) #5
      )
      # then
      flag = 1
    elseif (
        (k0==0) &&
        (lambda > -tol) && #1
        (k0 == 0 ? true : x0sort[k0] > theta + lambda - tol) && #2
        (theta + lambda >= x0sort[k0+1] - tol) && #3
        (x0sort[k1] >= theta - tol) && #4
        (k1 == n ? true : theta > x0sort[k1+1] - tol) #5
      )
      # then
      flag = 1
    elseif (
        (k1==n) &&
        (lambda > -tol) && #1
        (k0 == 0 ? true : x0sort[k0] > theta + lambda - tol) && #2
        (theta + lambda >= x0sort[k0+1] - tol) && #3
        # (x0sort[k1] >= theta - tol) && #4
        (x0sort[k1] - theta >= -tol) && #4
        (k1 == n ? true : theta > x0sort[k1+1] - tol) #5
      )
      # then
      flag = 1
    elseif (
        (lambda > -tol) && #1
        (k0 == 0 ? true : x0sort[k0] > theta + lambda - tol) && #2
        (theta + lambda >= x0sort[k0+1] - tol) && #3
        (x0sort[k1] >= theta - tol) && #4
        (k1 == n ? true : theta > x0sort[k1+1] - tol) #5
      )
      # then
      flag = 1
    end

    if flag == 1
      break
    elseif (flag == 0) && (k1 < n)
      k1 = k1 + 1
      sumk0p1k1 = sumk0p1k1 + x0sort[k1]
      if verb
        println("k1<n updates")
        println("x0sort[k1] = $(x0sort[k1])")
        println("sumk0p1k1 = $sumk0p1k1")
      end
    elseif (flag == 0) && (k0 > 0) && (k1 == n)
      k0 = k0 - 1
      k1 = k
      sumk0p1k += x0sort[k0+1]
      sumk0p1k1 = sumk0p1k
      sum1k0 -= x0sort[k0+1]
    end

    # tracking
    if hist
      history[nit+1] = (k0, k1)
    end
  end
end # time_pivot

  # primal
time_primal = @elapsed begin
  if x0prepop
    @simd for i in 1:k0
      @inbounds xbarsort[i] -= lambda
    end
  else
    @simd for i in 1:k0
      @inbounds xbarsort[i] = x0sort[i] - lambda
    end
  end
  @simd for i in k0+1:k1
    @inbounds xbarsort[i] = theta
  end
  if !x0prepop
    @simd for i in k1+1:n
      @inbounds xbarsort[i] = x0sort[i]
    end
  end
end # time_primal

  if verb
    println("k0 = $k0")
    println("k1 = $k1")
    println("θ = $theta")
    println("λ = $lambda")
    println("nit = $nit")
  end
  if hist
    return 1, (k0, k1), nit, (time_init_1+time_init_2, time_pivot, time_primal), hist, case, lambda
  else
    return 1, (k0, k1), nit, (time_init_1+time_init_2, time_pivot, time_primal), case, lambda
  end
end

#=================================================
GRB: detailed timing of GRBS and GRBU
  - x0sort: sorted input vector
  - r: scalar budget parameter
  - k: integer
=================================================#
include("../src/projection_gurobi.jl")