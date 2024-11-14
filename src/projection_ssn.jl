#=================================================
SSN method
  - b: sorted input vector (x0sort)
  - Πylambdapb: output (xbarsort)
  - lambda: top-k-sum functional
  - tau: scalar budget parameter (r)
  - nonneg: false=top-k-sum, true=vector-k-norm
  - verb: verbosity
  - prox_type: 1=stack, 2=pava; we observe that pava is much slower
  - tmp: additional workspace
=================================================#
function prox_stack!(x::AbstractVector{Tf}, lambda::Union{Tf,AbstractVector{Tf}}, nonneg::Bool) where {Tf<:Union{AbstractFloat,Rational}}
  """
  convenience
  """
  p = length(x)
  s = zeros(Tf, p)
  w = zeros(Tf, p)
  idx_i = zeros(Int64, p)
  idx_j = zeros(Int64, p)
  prox_stack!(x, s, w, idx_i, idx_j, lambda, nonneg)
end
function prox_stack!(x::AbstractVector{Tf}, s::AbstractVector{Tf}, w::AbstractVector{Tf}, idx_i::AbstractVector{Ti}, idx_j::AbstractVector{Ti}, lambda::Union{Tf,AbstractVector{Tf}}, nonneg::Bool) where {Tf<:Union{AbstractFloat,Rational},Ti<:Integer}
  """
  stack-based algorithm (Algorithm 4 in Bogdan et al. 2015)
  """
  p = length(x)
  # s = zeros(Tf, p)
  # w = zeros(Tf, p)
  d = zero(Tf)
  # idx_i = zeros(Int64, p)
  # idx_j = zeros(Int64, p)
  k = 1
  constlambda = length(lambda) == 1

  for i in 1:p
    idx_i[k] = i
    idx_j[k] = i
    s[k] = x[i] - (constlambda ? lambda : lambda[i])
    w[k] = s[k]

    while (k > 1) && (w[k - 1] <= w[k])
      k -= 1
      idx_j[k] = i
      s[k] += s[k + 1]
      w[k] = s[k] / (i - idx_i[k] + 1)
    end
    k += 1
  end

  if nonneg
    @simd for j in 1:k-1
      @inbounds d = max(w[j], 0.0)
      @simd for i in idx_i[j]:idx_j[j]
        @inbounds x[i] = d
      end
    end
  else
    @simd for j in 1:k-1
      @inbounds d = w[j]
      @simd for i in idx_i[j]:idx_j[j]
        @inbounds x[i] = d
      end
    end
  end
  nothing
end

function prox_pava!(y::AbstractVector{Tf}, lambda::Union{Tf,AbstractVector{Tf}}, nonneg::Bool) where {Tf<:Union{AbstractFloat,Rational}}
  """
  convenience
  """
  yc = zeros(Tf, n + 1)
  prox_pava!(y, yc, lambda, nonneg)
end
function prox_pava!(y::AbstractVector{Tf}, yc::AbstractVector{Tf}, lambda::Union{Tf,AbstractVector{Tf}}, nonneg::Bool) where {Tf<:Union{AbstractFloat,Rational}}
  """
  PAVA: modify y and yc
  """
  n = length(y)
  if lambda == zero(Tf)
    yc[2:end] .= cumsum(y)
  else
    yc[2:end] .= cumsum(y .- lambda)
  end
  known = 1
  tmp = zero(Tf)
  ip::Int64 = 1
  slope = -typemax(Tf)

  while known < n
    slope = -typemax(Tf)
    for i in known+1:n+1
      @inbounds tmp::Tf = (yc[i] - yc[known]) / (i - known)
      if tmp > slope
        slope = tmp
        ip = i
      end
    end
    @simd for i in known:ip-1
      @inbounds y[i] = (yc[ip] - yc[known]) / (ip - known)
    end
    known = ip
  end

  if nonneg
    y .= max.(y, 0.0)  # Ensure non-negativity
  end
  nothing
end

function Isupp!(Isupp::AbstractVector{Bool}, v::AbstractVector{Tf}, nonneg::Bool, tol::Tf=sqrt(eps)) where {Tf<:AbstractFloat}
  """
  Isupp(Π_C(d)) := {i : [B * Π_C(d)]_i == 0}
  """
  n = length(v)
  if nonneg; Isupp[n] = abs(v[n]) <= tol; end
  @simd for i in 1:n-1; @inbounds Isupp[i] = abs(v[i] - v[i+1]) <= tol; end
  nothing
end
function vHv!(v::AbstractVector{Tf}, Isupp::AbstractVector{Bool}) where {Tf<:AbstractFloat}
  """
  H is a projection matrix so v'Hv = (Hv)'Hv but here it is computed in single loop
  """
  n = length(v)
  if length(Isupp) != n
    throw(ArgumentError("The lengths of v and Isupp must be the same."))
  end

  # Initialize variables for computation
  current_sum = zero(Tf)
  current_count = 0
  start_idx = 1
  inner_product = zero(Tf)  # This will keep track of the inner product v' * Iv

  for i in 1:n
    # Add current value to the current sum
    current_sum += v[i]
    current_count += 1

    # Determine if the next index starts a new group or continues the current group
    if i == n || !Isupp[i]  # End of vector or start of a new group in the next step
      # Calculate the average for the current group
      group_average = current_sum / current_count

      # Calculate the part of the inner product for this group
      for fill_idx in start_idx:i
        inner_product += v[fill_idx] * group_average
      end

      # Reset for the next group
      start_idx = i + 1
      current_sum = zero(Tf)
      current_count = 0
    end
  end

  return inner_product  # Return the calculated inner product
end
function project_topksum_ssn!(
  Πylambdapb::AbstractVector{Tf},
  b::AbstractVector{Tf}, lambda::AbstractVector{Tf}, tau::Tf,
  nonneg::Bool=false, verb::Bool=false,
  prox_type::Ti=1, tmp=(zeros(Tf, n), zeros(Tf, n), zeros(Ti, n), zeros(Ti, n), zeros(Bool, n)), # prox_stack
) where {Tf<:Union{AbstractFloat,Rational},Ti<:Integer}
  """
  semismooth Newton method for solving the scalar equation ϕ'(y)=0, (3.7) in Li & Li 2022, where
    ϕ(y) := 0.5||Π_C(yλ + b)||^2 - yτ
    ϕ'(y) = ⟨Π_C(yλ + b), λ⟩ - τ
    C := {x : x_i ≥ x_i+1 ≥ 0} = {x : Bx ≥ 0}
    Π_C(v) := argmin 0.5||x-v||^2 : x ∈ C
    λ ∈ C
  """

  #
  # initialize
  #
  status = 0 # unsolved
  n = length(b)
  err = +Inf # error
  nit = 0 # newton iter
  lit = 0 # linesearch iter
  maxnit = 100
  maxlit = 6
  tol = 1e-10 # tolerance
  mu = 1e-3 # ls decrement
  delta = 0.9 # ls step scale
  
  y = 0.0 # iterate
  d = 0.0 # direction
  M = 0.0 # GJ
  phi_y = 1.0
  phiprime_y = 0.0
  phi_y_new = 0.0
  btb_2 = 0.5*(b'*b)
  Πylambdapb .= y .* lambda .+ b
  if prox_type == 1 # prox_stack!
    s_tmp = tmp[1]
    w_tmp = tmp[2]
    idx_i_tmp = tmp[3]
    idx_j_tmp = tmp[4]
    Isupp = tmp[5]
    @elapsed prox_stack!(Πylambdapb, s_tmp, w_tmp, idx_i_tmp, idx_j_tmp, zero(Tf), nonneg)
  elseif prox_type == 2 # prox_pava!
    yc_tmp = tmp[1]
    Isupp = tmp[2]
    prox_pava!(Πylambdapb, yc_tmp, zero(Tf), nonneg)
  else
    throw(ArgumentError("implemented prox solvers: prox_type:1 = prox_stack, prox_type:2 = prox_pava"))
  end
  
  #
  # iterate
  #
  while true
    
    # newton direction
    nit += 1
    Isupp!(Isupp, Πylambdapb, nonneg, tol)
    M = vHv!(lambda, Isupp)
    phiprime_y = Πylambdapb' * lambda - tau
    d = -phiprime_y / M

    err = abs(d)
    if err <= tol
      status = 1
      break
    elseif nit >= maxnit
      status = 2
      break
    end

    # linesearch
    phi_y = 0.5(Πylambdapb'*Πylambdapb) - y*tau - btb_2
    alpha = 1.0
    while true
      lit += 1
      y += alpha * d
      Πylambdapb .= y .* lambda .+ b
      if prox_type == 1
        prox_stack!(Πylambdapb, s_tmp, w_tmp, idx_i_tmp, idx_j_tmp, zero(Tf), nonneg)
      elseif prox_type == 2
        prox_pava!(Πylambdapb, yc_tmp, zero(Tf), nonneg)
      end
      phi_y_new = 0.5(Πylambdapb'*Πylambdapb) - y*tau - btb_2
      if phi_y_new <= phi_y + mu * alpha * phiprime_y || lit >= maxlit
        phi_y = phi_y_new
        break
      else
        alpha *= delta
      end
    end
  end
  return nit, lit, status
end