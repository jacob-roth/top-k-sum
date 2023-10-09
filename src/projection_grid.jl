function project_topksum_grid!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, maxt::Real=10_000,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
t_start = time()
  # inactive
  n = length(x0sort)
  if !active
    if !x0prepop
      @simd for i in 1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
    return 0, (-1, -1)
  end

  # initialize
  solved::Bool = false
  n = length(x0sort)
  s0::Tfr = sum(x0sort[1:k])
  lambda::Tfr = 0
  if k == n
    # xbarsort = x0sort - (s0 - r)/n
    lambda = (s0 - r) / k
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lambda
    end
    return 0, (0, n)
  elseif k == 1
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
    return 0, (0, findfirst(x0sort .< r))
  end

  #=
  preprocessing
  =#
  tol::Tfr = 0
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    tol = eps(Tfr)*x0sort[1]
  end
  n = length(x0sort)
  nit::Ti = 0
  
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
  
  # iterate
  while true
    if time() - t_start > maxt
      println("!TIMELIMIT!")
      return 0, (-1, -1)
    end
    nit += 1;
    kk0 = k - k0;
    k1k0 = k1 - k0;
    rho = k0 * k1k0 + kk0^2;
    
    sum1k0r_rho = (sum1k0 - r) / rho;
    sumk0p1k1_rho = sumk0p1k1 / rho;
    theta = k0 * sumk0p1k1_rho - kk0 * sum1k0r_rho;
    lambda = kk0 * sumk0p1k1_rho + k1k0 * sum1k0r_rho;
    
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
    elseif (flag == 0) && (k0 > 0) && (k1 == n)
      k0 = k0 - 1
      k1 = k
      sumk0p1k += x0sort[k0+1]
      sumk0p1k1 = sumk0p1k
      sum1k0 -= x0sort[k0+1]
    end
  end
  
  # primal
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
  return 1, (k0, k1)
end
