using FiniteDifferences
x = [1.,2.,3.,4.,5.]
c = [0.1, 0.5, 0.1, 0.2, 0.7]

g1(x) = sum(x.^2)/2
g2(x) = x[3]^3/3
g3(x) = c'*x
g4(x) = x[2]^3/3
g5(x) = x[1]^2/2
g6(x) = x[1]^3/3
G(x) = [g1(x); g2(x); g3(x); g4(x); g5(x); g6(x)]
JG(x) = jacobian(central_fdm(10, 1), G, x)[1]
nablag1(x) = grad(central_fdm(2, 1), g1, x)[1]
nablag2(x) = grad(central_fdm(2, 1), g2, x)[1]
nablag3(x) = grad(central_fdm(2, 1), g3, x)[1]
nablag4(x) = grad(central_fdm(2, 1), g4, x)[1]
nablag5(x) = grad(central_fdm(2, 1), g5, x)[1]
nablag6(x) = grad(central_fdm(2, 1), g6, x)[1]
nabla2g1(x) = jacobian(central_fdm(10, 1), nablag1, x)[1]
nabla2g2(x) = jacobian(central_fdm(10, 1), nablag2, x)[1]
nabla2g3(x) = jacobian(central_fdm(10, 1), nablag3, x)[1]
nabla2g4(x) = jacobian(central_fdm(10, 1), nablag4, x)[1]
nabla2g5(x) = jacobian(central_fdm(10, 1), nablag5, x)[1]
nabla2g6(x) = jacobian(central_fdm(10, 1), nablag6, x)[1]
proj(x) = begin
  xbar = similar(x)
  project_topksum_esgs2!(xbar, x, r, k, true)
  return xbar
end
Jproj(y) = jacobian(central_fdm(10, 1), proj, y)[1]

Gx = G(x)
Gxbar = similar(Gx)
r = 25.0
k = 3
varphi(x) = begin
  Gx = G(x)
  Gxbar = similar(Gx)
  project_topksum_esgs2!(Gxbar, Gx, r, k, true)
  return sum( (Gx .- Gxbar).^2 )/2
end

nablavarphi(x) = grad(central_fdm(2, 1), varphi, x)[1]
nabla2varphi(x) = jacobian(central_fdm(16, 1), nablavarphi, x)[1]
nablavarphi_manual(x) = begin
  Gx = G(x)
  Gxbar = similar(Gx)
  project_topksum_esgs2!(Gxbar, Gx, r, k, true)
  [
    nablag1(x)';
    nablag2(x)';
    nablag3(x)';
    nablag4(x)';
    nablag5(x)';
    nablag6(x)';
  ]' * (Gx - Gxbar)
end
nabla2varphi_manual(x) = begin
  Gx = G(x)
  Gxbar = similar(Gx)
  project_topksum_esgs2!(Gxbar, Gx, r, k, true)
  out = (
    + nabla2g1(x) * (Gx - Gxbar)[1]
    + nabla2g2(x) * (Gx - Gxbar)[2]
    + nabla2g3(x) * (Gx - Gxbar)[3]
    + nabla2g4(x) * (Gx - Gxbar)[4]
    + nabla2g5(x) * (Gx - Gxbar)[5]
    + nabla2g6(x) * (Gx - Gxbar)[6]
  )
  out += JG(x)' * (I - Jproj(Gx)) * JG(x)
  return out
end

hcat(nablavarphi(x), nablavarphi_manual(x))
norm(nablavarphi(x) - nablavarphi_manual(x), Inf)

nabla2varphi(x)
nabla2varphi_manual(x)
norm(nabla2varphi(x) .- nabla2varphi_manual(x), Inf)