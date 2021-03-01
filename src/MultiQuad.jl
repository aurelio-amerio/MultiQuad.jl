module MultiQuad

using QuadGK, Cuba, HCubature, Unitful, FastGaussQuadrature, LinearAlgebra, LRUCache

export quad, dblquad, tplquad


function _quadgk_wrapper_1D(arg::Function, x1, x2; kwargs...)
    return quadgk(arg, x1, x2; kwargs...)
end

function _cuba_wrapper_1D(arg::Function, x1, x2; method::Symbol, kwargs...)
    if method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
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

#order, method, optional

function get_weights_fastgaussquadrature(order::Int, method::Symbol, optional::Real=NaN)
    if method == :gausslegendre
        xw = gausslegendre
        return xw(order)

    elseif method == :gausshermite
        xw = gausshermite
        return xw(order)

    elseif method == :gausslaguerre
        xw = gausslaguerre
        if isnan(optional)
            return xw(order)
        else
            α = optional
            return xw(order, α)
        end

    elseif method == :gausschebyshev
        xw = gausschebyshev
        if isnan(optional)
            return xw(order)
        else
            kind = optional
            return xw(order, kind)
        end

    # elseif method == :gaussjacobi
    #     xw = gaussjacobi
    #     return xw(order)

    elseif method == :gaussradau
        xw = gaussradau
        return xw(order)

    elseif method == :gausslobatto
        xw = gausslobatto
        return xw(order)
    
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end

    
end

const lru_xw_fastgaussquadrature = LRU{
    Tuple{Int,Symbol,Real},
    Tuple{Vector{Float64},Vector{Float64}},
}(maxsize = Int(1e5))

function get_weights_fastgaussquadrature_cached(order::Int, method::Symbol, optional::Real=NaN)
    get!(lru_xw_fastgaussquadrature, (order, method, optional)) do
        get_weights_fastgaussquadrature(order, method, optional)
    end
end

function _fastgaussquadrature_wrapper_1D(arg::Function, x1, x2; method::Symbol, order::Int, kwargs...)
    
    if haskey(kwargs, :α)
        optional = kwargs[:α]
    elseif haskey(kwargs, :kind)
        optional = kwargs[:kind]
    else 
        optional = NaN
    end

    x, w = get_weights_fastgaussquadrature_cached(order, method, optional)

    b=x2
    a=x1

    b_a_2 = (b-a)/2
    a_plus_b_2 = (a+b)/2

    arg_gauss(x) = arg(b_a_2*x + a_plus_b_2)
    
    integral = b_a_2 * dot(w, arg_gauss.(x))
    return integral, NaN

end

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
function quad(arg::Function, x1, x2; method::Symbol = :quadgk, order::Int=1000, kwargs...)

    if method in [:suave, :vegas]
        return _cuba_wrapper_1D(arg, x1, x2; method=method)
    elseif method == :quadgk
        return _quadgk_wrapper_1D(arg, x1, x2; kwargs...)
    elseif method in [:gausslegendre, :gausshermite, :gausslaguerre, :gausschebyshev, :gaussjacobi, :gaussradau, :gausslobatto]
        return _fastgaussquadrature_wrapper_1D(arg, x1, x2; method=method, order=order, kwargs...)    
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end
end


function _cuba_wrapper_2D(arg::Function, x1, x2, y1::Function, y2::Function; method::Symbol, kwargs...)
    if method == :cuhre
        integrate = cuhre
    elseif method == :divonne
        integrate = divonne
    elseif method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end

    units = unit(arg(y1(x1), x1)) * unit(x1) * unit(y1(x1))

    arg1(a, x) = (y2(x) - y1(x)) * arg((y2(x) - y1(x)) * a + y1(x), x)


    arg2(a, b) = ustrip(units, (x2 - x1) * arg1(a, (x2 - x1) * b + x1))::Float64

    function integrand2(x, f)
        f[1] = arg2(x[1], x[2])
    end

    result, err = integrate(integrand2, 2, 1; kwargs...)

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end


function _hcubature_wrapper_2D(arg::Function, x1, x2, y1::Function, y2::Function; kwargs...)
    units = unit(arg(y1(x1), x1)) * unit(x1) * unit(y1(x1))

    arg1(a, x) = (y2(x) - y1(x)) * arg((y2(x) - y1(x)) * a + y1(x), x)


    arg2(a, b) = ustrip(units, (x2 - x1) * arg1(a, (x2 - x1) * b + x1))::Float64


    function integrand(arr)
        return arg2(arr[1], arr[2])
    end

    min_arr = [0, 0]
    max_arr = [1, 1]
    result, err = hcubature(integrand, min_arr, max_arr; kwargs...)
    

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

    if method in [:cuhre, :divonne, :suave, :vegas]
        return _cuba_wrapper_2D(arg, x1, x2, y1, y2; method=method)
    elseif method == :hcubature
        return _hcubature_wrapper_2D(arg, x1, x2, y1, y2; kwargs...)
    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end
end

function _cuba_wrapper_3D(arg::Function, x1, x2, y1::Function, y2::Function, z1::Function, z2::Function;  method::Symbol, kwargs...)
    if method == :cuhre
        integrate = cuhre
    elseif method == :divonne
        integrate = divonne
    elseif method == :suave
        integrate = suave
    elseif method == :vegas
        integrate = vegas
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

    
    function integrand2(x, f)
        f[1] = arg2(x[1], x[2], x[3])
    end

    result, err = integrate(integrand2, 3, 1; kwargs...)
    

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end

function _hcubature_wrapper_3D(arg::Function, x1, x2, y1::Function, y2::Function, z1::Function, z2::Function; kwargs...)

    units = unit(arg(z1(x1, y1(x1)), y1(x1), x1)) * unit(x1) * unit(y1(x1)) *
            unit(z1(y1(x1), x1))

    arg0(a, y, x) =
        (z2(x, y) - z1(x, y)) * arg((z2(x, y) - z1(x, y)) * a + z1(x, y), y, x)

    arg1(a, b, x) = (y2(x) - y1(x)) * arg0(a, (y2(x) - y1(x)) * b + y1(x), x)


    arg2(a, b, c) =
        ustrip(units, (x2 - x1) * arg1(a, b, (x2 - x1) * c + x1))::Float64

    
    function integrand(arr)
        return arg2(arr[1], arr[2], arr[3])
    end

    min_arr = [0, 0, 0]
    max_arr = [1, 1, 1]
    result, err = hcubature(integrand, min_arr, max_arr; kwargs...)
    
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
    method = :hcubature,
    kwargs...,
)

    if method in [:cuhre, :divonne, :suave, :vegas]
        return _cuba_wrapper_3D(arg, x1, x2, y1, y2, z1, z2;  method=method, kwargs...)

    elseif method == :hcubature
        return _hcubature_wrapper_3D(arg, x1, x2, y1, y2, z1, z2; kwargs...)

    else
        ex = ErrorException("Integration method $method is not supported!")
        throw(ex)
    end
end


end # module
