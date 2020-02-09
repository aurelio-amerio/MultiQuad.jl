module MultiQuad

using QuadGK, Cuba, HCubature, Unitful

export quad, dblquad, tplquad

@doc raw"""
    quad(arg::Function, x1, x2; method = :quadgk, kwargs...)

Performs the integral ``\int_{x1}^{x2}f(x)dx``

Available integration methods:
- `:suave`
- `:vegas`
- `:quadgk`

See [QuadGK](https://github.com/JuliaMath/QuadGK.jl) and [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

# Examples
```jldoctest
func(x) = x^2*exp(-x)
integral, error = quad(func, 0, 4)
```

`quad` can handle units of measurement through the [`Unitful`](https://github.com/PainterQubits/Unitful.jl) package:

```jldoctest
using Unitful

func(x) = x^2
integral, error = quad(func, 1u"m", 5u"m")
```
"""
function quad(arg::Function, x1, x2; method = :quadgk, kwargs...)

    if method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
    elseif method == :quadgk
        integrate = quadgk
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end

    if method == :quadgk
        return quadgk(arg, x1, x2; kwargs...)
    end

    units = unit(arg(x1)) * unit(x1)

    arg2(a) = ustrip(units, (x2 - x1) * arg((x2 - x1) * a + x1))::Float64

    function integrand(x, f)
        f[1] = arg2(x[1])
    end

    result, err = integrate(integrand, 1, 1; kwargs...)

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end

@doc raw"""
    dblquad(arg::Function, x1, x2, y1::Function, y2::Function; method = :cuhre, kwargs...)

Performs the integral ``\int_{x1}^{x2}\int_{y1(x)}^{y2(x)}f(y,x)dydx``

Available integration methods:
- `:cuhre`
- `:divonne`
- `:suave`
- `:vegas`
- `:hcubature`

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments fro the `:cuhre`, `:divonne`, `:suave` and `:vegas` methods.
See [HCubature](https://github.com/stevengj/HCubature.jl) for all the available keywords for the `:hcubature` method.

# Examples
```jldoctest
func(y,x) = sin(x)*y^2
integral, error = dblquad(func, 1, 2, x->0, x->x^2, rtol=1e-9)
```

`dblquad` can handle units of measurement through the [`Unitful`](https://github.com/PainterQubits/Unitful.jl) package:

```jldoctest
using Unitful

func(y,x) = x*y^2
integral, error = dblquad(func, 1u"m", 2u"m", x->0u"m^2", x->x^2, rtol=1e-9)
```
"""
function dblquad(
    arg::Function,
    x1,
    x2,
    y1::Function,
    y2::Function;
    method = :cuhre,
    kwargs...,
)

    if method == :cuhre
        integrate = cuhre
    elseif method == :divonne
        integrate = divonne
    elseif method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
    elseif method == :hcubature
        integrate = hcubature
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end

    units = unit(arg(y1(x1), x1)) * unit(x1) * unit(y1(x1))

    arg1(a, x) = (y2(x) - y1(x)) * arg((y2(x) - y1(x)) * a + y1(x), x)


    arg2(a, b) = ustrip(units, (x2 - x1) * arg1(a, (x2 - x1) * b + x1))::Float64

    if method == :hcubature
        function integrand(arr)
            return arg2(arr[1], arr[2])
        end

        min_arr = [0, 0]
        max_arr = [1, 1]
        result, err = integrate(integrand, min_arr, max_arr; kwargs...)
    else
        function integrand2(x, f)
            f[1] = arg2(x[1], x[2])
        end

        result, err = integrate(integrand2, 2, 1; kwargs...)
    end

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end

@doc raw"""
    tplquad(arg::Function, x1, x2, y1::Function, y2::Function, z1::Function, z2::Function; method = :cuhre, kwargs...)

Performs the integral ``\int_{x1}^{x2}\int_{y1(x)}^{y2(x)}\int_{z1(x,y)}^{z2(x,y)}f(z,y,x)dzdydx``

Available integration methods:
- `:cuhre`
- `:divonne`
- `:suave`
- `:vegas`
- `:hcubature`

See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments fro the `:cuhre`, `:divonne`, `:suave` and `:vegas` methods.
See [HCubature](https://github.com/stevengj/HCubature.jl) for all the available keywords for the `:hcubature` method.

# Examples
```jldoctest
func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0, 4, x->x, x->x^2, (x,y)->2, (x,y)->3*x)
```

`tplquad` can handle units of measurement through the [`Unitful`](https://github.com/PainterQubits/Unitful.jl) package:

```jldoctest
using Unitful

func(z,y,x) = sin(z)*y*x
integral, error = tplquad(func, 0u"m", 4u"m", x->0u"m^2", x->x^2, (x,y)->0, (x,y)->3)
```
"""
function tplquad(
    arg::Function,
    x1,
    x2,
    y1::Function,
    y2::Function,
    z1::Function,
    z2::Function;
    method = :cuhre,
    kwargs...,
)

    if method == :cuhre
        integrate = cuhre
    elseif method == :divonne
        integrate = divonne
    elseif method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
    elseif method == :hcubature
        integrate = hcubature
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end

    units = unit(arg(z1(x1, y1(x1)), y1(x1), x1)) * unit(x1) * unit(y1(x1)) *
            unit(z1(y1(x1), x1))

    arg0(a, y, x) =
        (z2(x, y) - z1(x, y)) * arg((z2(x, y) - z1(x, y)) * a + z1(x, y), y, x)

    arg1(a, b, x) = (y2(x) - y1(x)) * arg0(a, (y2(x) - y1(x)) * b + y1(x), x)


    arg2(a, b, c) =
        ustrip(units, (x2 - x1) * arg1(a, b, (x2 - x1) * c + x1))::Float64

    if method == :hcubature
        function integrand(arr)
            return arg2(arr[1], arr[2], arr[3])
        end

        min_arr = [0, 0, 0]
        max_arr = [1, 1, 1]
        result, err = integrate(integrand, min_arr, max_arr; kwargs...)
    else
        function integrand2(x, f)
            f[1] = arg2(x[1], x[2], x[3])
        end

        result, err = integrate(integrand2, 3, 1; kwargs...)
    end

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end


end # module
