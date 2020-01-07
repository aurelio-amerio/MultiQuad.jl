module MultiQuad

@doc raw"""
    dblquad(arg::Function, x1, x2, y1::Function, y2::Function; method = :cuhre, kwargs...)
Performs the integral ``\int_{x1}^{x2}\int_{y1(x)}^{y2(x)}f(y,x)dydx``
See [Cuba.jl](https://giordano.github.io/Cuba.jl/stable/) for all the available keyword arguments.

# Examples
```jldoctest
func(y,x) = sin(x)*y^2
dblquad(func, 1, 2, y->0, y->y^2, rtol=1e-9)
```

`dblquad` can handle units of measurement through the [`Unitful`](https://github.com/PainterQubits/Unitful.jl) package:

```jldoctest
using Unitful

func(y,x) = sin(x)*y^2
dblquad_cuhre(func, 1u"m", 2u"m", y->0, y->y^2, rtol=1e-9)
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
    else
        ex = ErrorException("The selected integration method $method is not supported!")
        throw(ex)
    end

    units = unit(arg(y1(x1), x1)) * unit(x1) * unit(y1(x1))

    arg1(a, x) = (y2(x) - y1(x)) * arg((y2(x) - y1(x)) * a + y1(x), x)


    arg2(a, b) =
        ustrip(units, (x2 - x1) * arg1(a, (x2 - x1) * b + x1))::Float64

    function integrand(x, f)
        f[1] = arg2(x[1], x[2])
    end

    result, err = integrate(integrand, 2, 1; kwargs...)

    if units == Unitful.NoUnits
        return result[1], err[1]
    else
        return (result[1], err[1]) .* units
    end
end

# function tplquad(
#     arg::Function,
#     x1,
#     x2,
#     y1::Function,
#     y2::Function.
#     z1::Function,
#     z2::Function;
#     method = :cuhre,
#     kwargs...,
# )
#
#     if method == :cuhre
#         integrate = cuhre
#     elseif method == :divonne
#         integrate = divonne
#     elseif method == :suave
#         integrate = suave
#     elseif method == :vegas
#         integrate = vegas
#     else
#         ex = ErrorException("The selected integration method $method is not supported!")
#         throw(ex)
#     end
#
#     units = unit(arg(z1(y1(x)), y1(x1), x1)) * unit(x1) * unit(y1(x1)) * unit(z1(y1(x1),x1))
#
#     arg0(a, y, x) = (z2(y2(x), x) - z1(x))
#
#     arg1(a, x) = (y2(x) - y1(x)) * arg((y2(x) - y1(x)) * a + y1(x), x)
#
#
#     arg2(a, b) =
#         ustrip(units, (x2 - x1) * arg1(a, (x2 - x1) * b + x1))::Float64
#
#     function integrand(x, f)
#         f[1] = arg2(x[1], x[2])
#     end
#
#     result, err = integrate(integrand, 2, 1; kwargs...)
#
#     if units == Unitful.NoUnits
#         return result[1], err[1]
#     else
#         return (result[1], err[1]) .* units
#     end
# end


end # module
