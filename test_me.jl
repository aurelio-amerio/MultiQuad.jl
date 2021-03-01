using Pkg
Pkg.activate("./")
Pkg.instantiate()
# using Revise
using MultiQuad
using Unitful

#%%
f(x) = x^4


@test quad(f, -1, 1, method=:gausslegendre, order=100000)[1] ≈ 2/5

@test quad(f, -1, 1, method=:gausshermite, order=10000)[1] ≈ 3(√π)/4

@test quad(f, -1, 1, method=:gausslaguerre, order=10000)[1] ≈ 24
@test quad(f, -1, 1, method=:gausslaguerre, order=3, α=1.0)[1] ≈ 120

@test quad(f, -1, 1, method=:gausschebyshev, order=3)[1] ≈ 3π/8
@test quad(f, -1, 1, method=:gausschebyshev, order=3, kind=1)[1] ≈ 3π/8
@test quad(f, -1, 1, method=:gausschebyshev, order=3, kind=2)[1] ≈ π/16
@test quad(f, -1, 1, method=:gausschebyshev, order=3, kind=3)[1] ≈ 3π/8
@test quad(f, -1, 1, method=:gausschebyshev, order=3, kind=4)[1] ≈ 3π/8

@test quad(f, -1, 1, method=:gaussradau, order=3)[1] ≈ 2/5
@test quad(f, -1, 1, method=:gausslobatto, order=4)[1] ≈ 2/5




#%%
