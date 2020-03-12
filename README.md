[![Build Status](https://travis-ci.com/aurelio-amerio/MultiQuad.jl.svg?branch=master)](https://travis-ci.com/aurelio-amerio/MultiQuad.jl)
[![Coverage Status](https://coveralls.io/repos/github/aurelio-amerio/MultiQuad.jl/badge.svg?branch=master)](https://coveralls.io/github/aurelio-amerio/MultiQuad.jl?branch=master)
[![codecov.io](https://codecov.io/github/aurelio-amerio/MultiQuad.jl/coverage.svg?branch=master)](https://codecov.io/github/aurelio-amerio/MultiQuad.jl?branch=master)

# MultiQuad.jl
**MultiQuad.jl** is a convenient interface to perform 1D, 2D and 3D numerical integrations.
It uses [QuadGK](https://github.com/JuliaMath/QuadGK.jl), [Cuba](https://github.com/giordano/Cuba.jl) and [HCubature](https://github.com/stevengj/HCubature.jl) as back-ends.

# Usage

## `quad`

```julia
quad(arg::Function, x1, x2; method = :quadgk, kwargs...)
```

It is possible to use `quad` to perform 1D integrals of the following kind:

<img src="https://quicklatex.com/cache3/1c/ql_cb17b9863c232f8f8b8219f642d7771c_l3.png" title="\int_{x1}^{x2}f(x)dx" />

The supported integration method are:

- `:quadgk`
- `:vegas`
- `:suave`

It is suggested to use `:quadgk` (the default option).

See [QuadGK](https://github.com/JuliaMath/QuadGK.jl) and [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://quicklatex.com/cache3/d0/ql_79904108afd5b30c1ba38a483e6d93d0_l3.png" title="\int_0^4 x^2e^{-x}dx" />

```julia
using MultiQuad

func(x) = x^2*exp(-x)
quad(func, 0, 4)
```

### Example 2

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://quicklatex.com/cache3/67/ql_3183fbe4d2c967dc939a860aff773b67_l3.png" title="\int_{0 m}^{5 m} x^2e^{-x}dx" />

```julia
using MultiQuad
using Unitful

func(x) = x^2
quad(func, 1u"m", 5u"m")
```

## `dblquad`

```julia
dblquad(arg::Function, x1, x2, y1::Function, y2::Function; method = :cuhre, kwargs...)
```

It is possible to use `quad` to perform 2D integrals of the following kind:

<img src="https://quicklatex.com/cache3/98/ql_4ff96f9007a0dab7c9e474074729c998_l3.png" title="\int_{x1}^{x2}\int_{y1(x)}^{y2(x)}f(y,x)dydx" />

The supported integration method are:

- `:cuhre` (default)
- `:vegas`
- `:suave`
- `:divonne`

It is suggested to use `:cuhre` (the default option).

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://quicklatex.com/cache3/14/ql_f5c7d790447fe9b2b85baab592aecb14_l3.png" title="\int_1^2\int_{0}^{x^2}\sin(x) y^2 dxdy" />

```julia
using MultiQuad

func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1, 2, x->0, x->x^2, rtol=1e-9)
```

### Example 2 

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://quicklatex.com/cache3/2f/ql_5a32b5f5f6cd8b24ca4094702e27562f_l3.png" title="\int_{1m}^{2m}\int_{0m^2}^{x^2}x y^2 dxdy" />

```julia
using MultiQuad
using Unitful

func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1u"m", 2u"m", x->0u"m^2", x->x^2, rtol=1e-9)
```

## `tplquad`

```julia
tplquad(arg::Function, x1, x2, y1::Function, y2::Function, z1::Function, z2::Function; method = :cuhre, kwargs...)
```

It is possible to use `quad` to perform 3D integrals of the following kind:

<img src="https://quicklatex.com/cache3/a2/ql_99483faa353fa0686da01c885323eba2_l3.png" title="\int_{x1}^{x2}\int_{y1(x)}^{y2(x)}\int_{z1(x,y)}^{z2(x,y)}f(z,y,x)dzdydx" />

The supported integration method are:

- `:cuhre` (default)
- `:vegas`
- `:suave`
- `:divonne`

It is suggested to use `:cuhre` (the default option)

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://quicklatex.com/cache3/3f/ql_1405adbfe074a85b1e98ddc359b7553f_l3.png" title="\int_{0}^{4}\int_{x}^{x^2}\int_{2}^{3x}\sin(z) \, x \, y \, dzdydx" />

```julia
using MultiQuad

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0, 4, x->x, x->x^2, (x,y)->2, (x,y)->3*x)
```

### Example 2

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://quicklatex.com/cache3/a3/ql_246a3dfa043247febb3aeca76d30efa3_l3.png" title="\int_{0m}^{4m}\int_{0m^2}^{x^2}\int_{0}^{3}\sin(z) \, x \, y \, dzdydx" />

```julia
using MultiQuad
using Unitful

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0u"m", 4u"m", x->0u"m^2", x->x^2, (x,y)->0, (x,y)->3)
```

