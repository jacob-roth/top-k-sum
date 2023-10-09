"""
(slow) max-k-sum calc
"""
function maxksum(x0::Vector, k::Integer)
  return sum(x0[sortperm(x0, rev=true)][1:k])
end

"""
compute indices `k0` and `k1` for a sorted vector with given `k`
"""
function get_k0k1(xsort::AbstractVector{T}, k::Int64, tol=eps()) where T<:Union{Rational,AbstractFloat}
  n = length(xsort)
  rtol = tol*max(1.0, maximum(xsort))
  # @assert(issorted(xsort, rev=true, lt=(x,y)->(x<y-rtol)))
  #! @assert(issorted(xsort, rev=true, lt=<))
  # k0 = findlast(xsort .> xsort[k] + rtol) # good
  k0 = findlast(xsort .> xsort[k])
  # k0 = findlast(xsort .> xsort[k] - rtol) # bad
  k0 = isnothing(k0) ? 0 : k0
  # k1 = findfirst(xsort .< xsort[k] - rtol) # good
  k1 = findfirst(xsort .< xsort[k])
  # k1 = findfirst(xsort .< xsort[k] + rtol) # bad
  k1 = isnothing(k1) ? n : k1 - 1
  return k0, k1
end
function get_k0k1!(xsort::AbstractVector{T}, k::Int64, tol=eps()) where T<:Union{Rational,AbstractFloat}
  n = length(xsort)
  rtol = tol*max(1.0, maximum(xsort))
  # @assert(issorted(xsort, rev=true, lt=(x,y)->(x<y-rtol))) # approximately sorted; #! removed for speed
  # @assert(issorted(xsort, rev=true, lt=<)) # actually sorted; #! removed for speed
  k0 = 0
  k1 = n
  xsortk = xsort[k]
  rhs = xsortk + rtol
  @inbounds for i in eachindex(xsort)
    if xsort[i] > rhs
      # new: if-elseif-end clause
      # safeguard against numerical error in sorting
      if i == 1 && xsort[i+1] <= rhs # safeguard i==1
        k0 = i
      elseif i > 1 && i < n && xsort[i-1] > rhs && xsort[i+1] <= rhs # safeguard i∈[2,n-1]
        k0 = i
      elseif i == n && xsort[i-1] > rhs # safeguard i==n
        k0 = i
      end
      # previously: just had below line
      # k0 = i
    else
      # if (xsort[i] <= xsortk + rtol) && (xsort[i] >= xsortk - rtol)
      if isapprox(xsort[i], xsortk, rtol=tol)
        continue
      elseif k1 == n # xsort[i] < xsortk
        k1 = i-1
      end
    end
  end
  return k0, k1
end

# """
# A matrix in conjugate function of maxksum indicator
# δ⁺(z) = max_y { ⟨y , Amat⁻ᵀ zꜜ⟩ : y <= [r; 0_{n-1}] }
# """
# function Amat!(A::Mf, n::Ti, k::Ti) where {Ti<:Integer, Tf<:AbstractFloat, Mf<:AbstractMatrix{Tf}}
#   A[1,:] = [ones(Tf, k)' zeros(Tf, n-k)']
#   A[2:end,:] .= -diagm(0=>ones(Tf, n), 1=>-ones(Tf, n))[1:n-1,1:n]
#   A
# end
# function Amat(n::Ti, k::Ti, Tf::DataType=Float64) where {Ti<:Integer}
#   A = zeros(Tf, n, n)
#   Amat!(A, n, k)
# end
# function kAmatinv!(kAinv::Mf, n::Ti, k::Ti) where {Ti<:Integer, Tf<:Real, Mf<:AbstractMatrix{Tf}}
#   # Ainv a float matrix (https://stackoverflow.com/questions/45373679/why-is-it-faster-to-perform-float-by-float-matrix-multiplication-compared-to-int)
#   for i in 1:n
#     for j in 1:i
#       kAinv[i,j] = max(one(Tf), min(j-1, k))
#     end
#     for j in i+1:n
#       kAinv[i,j] = min(j-k-1, zero(Tf))
#     end
#   end
#   kAinv
# end
# function Amatinv(n::Ti, k::Ti, Tf::DataType=Float64) where {Ti<:Integer}
#   Ainv = zeros(Tf, n, n)
#   kAmatinv!(Ainv, n, k)
#   Ainv .*= 1/k
# end
# struct Amatinv_action{Ti<:Integer} <: Function
#   n::Ti
#   k::Ti
# end
# function (B::Amatinv_action)(Bv::AbstractVector{Tf}, v::AbstractVector{Tf}) where Tf<:AbstractFloat
#   j_hi::Int64 = 0
#   @inbounds for i in 1:B.n
#     Bv[i] = 0.0
#     @inbounds for j in 1:i
#       Bv[i] += max(one(Tf), min(j-1, B.k)) * v[j] / B.k
#     end
#     j_hi = max(i+1, B.k)
#     @inbounds for j in i+1:j_hi
#       Bv[i] += min(j-B.k-1, zero(Tf)) * v[j] / B.k
#     end
#   end
#   nothing
# end
# struct Atmatinv_action{Ti<:Integer} <: Function
#   n::Ti
#   k::Ti
# end
# function (Bt::Atmatinv_action)(Btv::AbstractVector{Tf}, v::AbstractVector{Tf}) where Tf<:AbstractFloat
#   reverse!(v)
#   cumsum!(Btv, v)
#   reverse!(Btv)
#   reverse!(v)
#   Btv[1] *= 1/Bt.k
#   s = 0.0
#   @simd for i in 2:Bt.k
#     @inbounds Btv[i] *= (i-1)/Bt.k
#     # @inbounds Btv[i] -= sum(view(v, 1:i-1)) * (1-(i-1)/Bt.k)
#     s += v[i-1]
#     @inbounds Btv[i] -= s * (1-(i-1)/Bt.k)
#   end
#   nothing
# end

function had!(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T<:Number}
  """hadamard"""
    m,n = size(A)
    @assert (m,n) == size(B)
    for j in 1:n
       for i in 1:m
         @inbounds A[i,j] *= B[i,j]
       end
    end
    return A
end
function had!(A::AbstractVector{T}, B::AbstractVector{T}) where {T<:Number}
  """hadamard"""
    n = length(A)
    @assert n == length(B)
    @simd for i in 1:n
      @inbounds A[i] *= B[i]
    end
    return A
end

function get_thetalambda(xsort::Vector{T}, x0sort::Vector{T}, r::T, k::Int64) where T<:Union{Rational,AbstractFloat}
  @assert(issorted(xsort, rev=true))
  @assert(issorted(x0sort, rev=true))
  k0, k1 = get_k0k1(xsort, k)
  rho = k0 * (k1 - k0) + (k - k0)^2
  theta = (k0 * sum(x0sort[k0+1:k1]) - (k - k0) * (sum(x0sort[1:k0]) - r)) / rho
  lambda = ((k - k0) * sum(x0sort[k0+1:k1]) + (k1 - k0) * (sum(x0sort[1:k0]) - r)) / rho
  return rho, theta, lambda
end
function get_mu(xsort::Vector{T}, x0sort::Vector{T}, r::T, k::Int64) where T<:Union{Rational,AbstractFloat}
  @assert(issorted(xsort, rev=true))
  @assert(issorted(x0sort, rev=true))
  if sum(x0sort[1:k]) <= r
    mu = zeros(length(x0sort))
  else
    rho, theta, lambda = get_thetalambda(xsort, x0sort, r, k)
    mu = (x0sort - xsort) ./ lambda
  end
  return mu
end

#
# lcp helpers
#

function minv(idx, jdx, alpha)
  nalpha = length(alpha)
  return (((nalpha+1) - max(idx,jdx)) * min(idx,jdx)) // (nalpha+1)
end
function get_q(x0sort)
  return [x0sort[i] - x0sort[i+1] for i in 1:length(x0sort)-1]
end
function get_iM(i,j,absalpha)
  return (((absalpha+1) - max(i,j)) * min(i,j)) // (absalpha+1)
end
function get_M(alpha)
  absalpha = length(alpha)
  return diagm(0=>2//1 .* ones(Int64, absalpha), 1=>-1//1 .* ones(Int64, absalpha-1), -1=>-1//1 .* ones(Int64, absalpha-1))
end
function get_B(alpha, n)
  return Int.(diagm(0=>ones(n), 1=>-ones(n))[1:n-1,1:n]) .// 1
end

#
# grid helpers
#

function lambda(k0, k1, k, r, x0sort)
  num = (k-k0) * sum(x0sort[k0+1:k1]) + (k1-k0) * (sum(x0sort[1:k0]) - r)
  rho = k0*(k1-k0) + (k-k0)^2
  return num / rho
end
function theta(k0, k1, k, r, x0sort)
  num = k0 * sum(x0sort[k0+1:k1]) - (k - k0) * (sum(x0sort[1:k0]) - r)
  rho = k0*(k1-k0) + (k-k0)^2
  return num / rho
end