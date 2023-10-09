function ypBtz!(
  y::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr},
  a::Ti, b::Ti, k_alpha::Ti, lam::Tfr
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
  # initialize
  m = b - a + 1
  c = zero(Tfr)
  ai = zero(Ti)
  zer0 = zero(Tfr)

  # alpha[1]
  for i in 1:m
    ai = b-i+1
    @inbounds c += i * (
      -(x0sort[ai] - x0sort[ai+1]) + (
        m-i+1 == k_alpha ? +lam : zer0
      )
    )
  end
  c /= (m+1)
  y[a] += c

  # alpha[2],...,alpha[m]+1
  for i in 2:m+1
    ai = a+i-2
    @inbounds c -= (
      -(x0sort[ai] - x0sort[ai+1]) + (
        i-1 == k_alpha ? +lam : zer0
      )
    )
    y[a+i-1] += c
  end
  nothing
end

function project_maxksum_plcp!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool=true, x0prepop::Bool=false
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
  t::Ti = 0
  s0::Tfr = sum(view(x0sort, 1:k))
  qk::Tfr = (k == n ? typemax(Tfr) : x0sort[k] - x0sort[k+1])
  lam::Tfr = 0
  
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
  else
    # iteration t = 0
    # evaluate mks @ x(lam_1): lam_1=q_k & z=0
    lam = qk
    m = s0 - k * lam
    if m <= r
      lam = (s0 - r) / k
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
    if solved
      return 0, (k-1, k)
    end
  end

  # iteration t = 1
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

  # iterate
  while true
    # get breakpoint
    lam_a = (a > 1 ? ((x0sort[a-1] - x0sort[a]) - za) / minv_ak : pInf) # (q[a-1] - za) / minv_ka
    if b+1 < n
      lam_b = ((x0sort[b+1] - x0sort[b+2]) - zb) / minv_bk
    else
      lam_b = pInf
    end
    lam = min(lam_a, lam_b)
    
    # check max-k-sum
    if lam != pInf
      m = s0 - k * lam + zk + lam * minv_kk
    else
      m = nInf
    end
    if m <= r# + tol
      lam = (s0 - r + zk) / (k - minv_kk)
      solved = true
      break
    end
    
    # update solution
    if lam_a <= lam_b
      # s = a-1
      za = (za - (x0sort[a-1] - x0sort[a])) / sigma # q[a-1] = x0sort[a-1] - x0sort[a]
      zk = zk + za * minv_ak
      zb = zb + za * minv_ab
      a -= 1
      k_alpha += 1
      b_alpha += 1
    else
      # s = b+1
      zb = (zb - (x0sort[b+1] - x0sort[b+2])) / sigma # q[b+1] = x0sort[b+1] - x0sort[b+2]
      za = za + zb * minv_ab
      zk = zk + zb * minv_bk
      b += 1
      b_alpha += 1
    end

    # increment basis length
    t += 1

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
  end

  if solved
    if !x0prepop
      xbarsort .= x0sort
    end
    for i in 1:k
      xbarsort[i] -= lam
    end
    ypBtz!(xbarsort, x0sort, a, b, k_alpha, lam)
  end

  # return
  return 1, (max(a-1, 0), min(b+1, n))
end