function get_isotonic(n::Ti, k::Ti, DT::DataType=Rational{Ti}) where Ti<:Integer
  """eq. 7: Bx <= b"""
  [
    [ones(DT, k)' zeros(DT, n-k)'];
    -diagm(0=>ones(DT, n), 1=>-ones(DT, n))[1:n-1,1:n] # -D
  ]
end
function get_isotonic_moreau(n::Ti, k::Ti, DT::DataType=Rational{Ti}) where Ti<:Integer
  Binvt = zeros(DT, n, n)
  for i in 1:n
    if i == 1
      Binvt[1,:] .= +1
    elseif i <= k
      Binvt[i,1:i-1] .= -(k-i+1)
      Binvt[i,i:end] .= +(i-1)
    else
      Binvt[i,i:end] .= + k
    end
  end
  return Binvt.//k
end
function get_unsortedtopk(x0::Vector, k::Integer)
  """eq. 5: Bx <= b"""
  n = length(x0)
  x0sort = sort(x0, rev=true)
  x0_k_ = x0sort[k]
  _k_ = findfirst(x0 .== x0_k_)
  kappa = findall(x0 .>= x0_k_)[1:k]
  Bprime = zeros(n, n)
  for i in 1:n
    if i in kappa # Bprime_I1
      Bprime[i,i] = -1
      Bprime[i,_k_] = +1
    else # Bprime_I2
      Bprime[i,i] = +1
      Bprime[i,_k_] = -1
    end
  end
  kappa_minus_k_ = setdiff(kappa, _k_)
  kappac = setdiff(collect(1:n), kappa)
  indkappa = zeros(n)
  indkappa[kappa] .= +1.0
  B = [indkappa'; Bprime[kappa_minus_k_,:]; Bprime[kappac,:]]
  return B
end

function get_psortedtopk(x0::Vector, k::Integer)
  """eq. 5: Bx <= b"""
  @assert(issorted(x0[1:k]))
  get_unsortedtopk(x0, k)  
end
