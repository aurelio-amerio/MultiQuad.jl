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

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{x1}^{x2}f(x)dx" title="\int_{x1}^{x2}f(x)dx" />

The supported adaptive integration methods are:

- `:quadgk`
- `:vegas`
- `:suave`

There are several fixed order quadrature methods available based on [`FastGaussQuadrature`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)

- `:gausslegendre` 
- `:gausshermite` 
- `:gausslaguerre` 
- `:gausschebyshev` 
- `:gaussradau` 
- `:gausslobatto`

If a specific integration routine is not needed, it is suggested to use `:quadgk` (the default option).
For highly oscillating integrals, it is advised to use `:gausslegendre`with a high order (~10000).

Please note that `:gausshermite` and `:gausslaguerre` are used to specific kind of integrals with infinite bounds.

For more informations on these integration techniques, please follow the official documentation [here](https://juliaapproximation.github.io/FastGaussQuadrature.jl/stable/gaussquadrature/).

See [QuadGK](https://github.com/JuliaMath/QuadGK.jl) and [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{0}^{4}&space;x^2e^{-x}dx" title="\int_{0}^{4} x^2e^{-x}dx" />

```julia
using MultiQuad

func(x) = x^2*exp(-x)
quad(func, 0, 4) # equivalent to quad(func, 0, 4, method=:quadgk)
quad(func, 0, 4, method=:gausslegendre, order=1000)

# for certain kinds of integrals with infinite bounds, it may be possible to use a specific integration routine
func(x) = x^2*exp(-x)
quad(func, 0, Inf, method=:quadgk)[1]
# gausslaguerre computes the integral of f(x)*exp(-x) from 0 to infinity
quad(x -> x^4, method=:gausslaguerre, order=10000)[1] 
```

### Example 2

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{0&space;m}^{5&space;m}&space;x^2e^{-x}dx" title="\int_{0 m}^{5 m} x^2e^{-x}dx" />

```julia
using MultiQuad
using Unitful

func(x) = x^2
quad(func, 1u"m", 5u"m")
```

## `dblquad`

```julia
dblquad(arg::Function, x1, x2, y1::Function, y2::Function; method = :cuhre, kwargs...)
dblquad(arg::Function, x1, x2, y1, y2; method = :cuhre, kwargs...)
```

It is possible to use `dblquad` to perform 2D integrals of the following kind:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{x1}^{x2}dx\int_{y1(x)}^{y2(x)}dyf(y,x)" title="\int_{x1}^{x2}dx\int_{y1(x)}^{y2(x)}dyf(y,x)" />

The supported integration method are:

- `hcubature` (default)
- `:cuhre` 
- `:vegas`
- `:suave`
- `:divonne`

It is suggested to use `:hcubature` (the default option).

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) and [HCubature](https://github.com/stevengj/HCubature.jl) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_1^2&space;dx&space;\int_{0}^{x^2}dy&space;\sin(x)&space;y^2" title="\int_1^2 dx \int_{0}^{x^2}dy \sin(x) y^2" />

```julia
using MultiQuad

func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1, 2, x->0, x->x^2, rtol=1e-9)
```

### Example 2 

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{1m}^{2m}dy\int_{0m^2}^{x^2}&space;dx&space;\,&space;x&space;y^2" title="\int_{1m}^{2m}dy\int_{0m^2}^{x^2} dx \, x y^2" />

```julia
using MultiQuad
using Unitful

func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1u"m", 2u"m", x->0u"m^2", x->x^2, rtol=1e-9)
```

### Example 3

Compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_1^2&space;dx&space;\int_{0}^{4}dy&space;\sin(x)&space;y^2" title="\int_1^2 dx \int_{0}^{4}dy \sin(x) y^2" />

```julia
using MultiQuad

func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1, 2, 0, 4, rtol=1e-9)
```

## `tplquad`

```julia
tplquad(arg::Function, x1, x2, y1::Function, y2::Function, z1::Function, z2::Function; method = :cuhre, kwargs...)
```

It is possible to use `quad` to perform 3D integrals of the following kind:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{x1}^{x2}dx\int_{y1(x)}^{y2(x)}dy\int_{z1(x,y)}^{z2(x,y)}dz&space;\,&space;f(z,y,x)" title="\int_{x1}^{x2}dx\int_{y1(x)}^{y2(x)}dy\int_{z1(x,y)}^{z2(x,y)}dz \, f(z,y,x)" />

The supported integration method are:

- `:cuhre` (default)
- `:vegas`
- `:suave`
- `:divonne`

It is suggested to use `:cuhre` (the default option)

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

### Example 1

Compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{0}^{4}dx\int_{x}^{x^2}dy\int_{2}^{3x}dz\sin(z)&space;\,&space;x&space;\,&space;y" title="\int_{0}^{4}dx\int_{x}^{x^2}dy\int_{2}^{3x}dz\sin(z) \, x \, y" />

```julia
using MultiQuad

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0, 4, x->x, x->x^2, (x,y)->2, (x,y)->3*x)
```

### Example 2

It is possible to compute integrals with unit of measurement using `Unitful`. 

For example, let's compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{0m}^{4m}dx\int_{0m^2}^{x^2}dy\int_{0}^{3}dz\sin(z)&space;\,&space;x&space;\,&space;y" title="\int_{0m}^{4m}dx\int_{0m^2}^{x^2}dy\int_{0}^{3}dz\sin(z) \, x \, y" />

```julia
using MultiQuad
using Unitful

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0u"m", 4u"m", x->0u"m^2", x->x^2, (x,y)->0, (x,y)->3)
```

### Example 3

Compute:

<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\int_{0}^{4}dx\int_{x}^{2}dy\int_{2}^{3}dz\sin(z)&space;\,&space;x&space;\,&space;y" title="\int_{0}^{4}dx\int_{1}^{2}dy\int_{2}^{3}dz\sin(z) \, x \, y" />

```julia
using MultiQuad

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0, 4, 1, 2, 2, 3)
```