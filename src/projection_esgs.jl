function project_maxksum_esgs!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
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
  s0::Tfr = sum(view(x0sort, 1:k))
  lam::Tfr = 0
  case::Ti = 0

  if k == n
    # xbarsort = x0sort - (s0 - r)/n
    lam = (s0 - r) / k
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lam
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

  # preprocessing
  tol = zero(Tfr)
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    tol = eps(Tfr)*x0sort[1]
  end
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
  nit::Ti = 0
  sum1k0::Tfr = s0 - x0sort[k]
  sumk0p1k1::Tfr = x0sort[k]
  sum1k0r = zero(Tfr)
  rho = zero(Tfr)
  solved = false
  
  # iterate
  while true

    # counter
    nit += 1

    # rho, theta, and lambda
    kk0 = k - k0
    k1k0 = k1 - k0
    rho = k0 * k1k0 + kk0^2
    sum1k0r = sum1k0 - r
    the = (k0 * sumk0p1k1 - kk0 * sum1k0r) / rho
    if k0 > 0
      theplam = (k * the + sum1k0r) / k0
    else
      theplam = (k * sumk0p1k1 + (k1-k) * sum1k0r) / rho
    end

    # kkt conditions
    kkt2 = (k0 == 0 ? true : x0sort[k0] > theplam-tol)
    kkt5 = (k1 == n ? true : the > x0sort[k1+1]-tol)

    # move k0
    if (kkt2 && kkt5) || (k0 == 0 && k1 == n)
      solved = true
      break
    elseif kkt2 && (k1 < n)
      k1 += 1
      sumk0p1k1 += x0sort[k1]
    elseif !kkt2 && (k0 > 0)
      sum1k0 -= x0sort[k0]
      sumk0p1k1 += (k0 > 0 ? x0sort[k0] : 0)
      k0 -= 1
    end
  end

  # stop
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
  if solved
    return 1, (k0, k1)
  end
end
