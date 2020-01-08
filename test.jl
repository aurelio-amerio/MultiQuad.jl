using MultiQuad, Test, Unitful
rtol=1e-10
xmin = 0
xmax = 4
func1(x) = x^2 * exp(-x)
res = 2 - 26/exp(4)
@test abs(quad(func1, xmin, xmax, rtol = rtol, method=:vegas)[1] - res) < 1e-5
@test abs(quad(func1, xmin, xmax, atol=1e-14, method=:suave)[1] - res) < 1e-3

#%%
rtol = 1e-10

xmin = 0
xmax = 4
ymin(x) = x^2
ymax(x) = x^3

res=241664/5
func1(y, x) = y * x^3
@test abs(dblquad(func1, xmin, xmax, ymin, ymax, rtol = rtol)[1] - res) < rtol * 1e2
@test abs(dblquad(func1, xmin, xmax, ymin, ymax, rtol = rtol, method=:vegas)[1] - res) < rtol * 1e5
@test abs(dblquad(func1, xmin, xmax, ymin, ymax, rtol = rtol, method=:suave)[1] - res) < rtol * 1e5
@test abs(dblquad(func1, xmin, xmax, ymin, ymax, rtol = rtol, method=:divonne)[1] - res) < rtol * 1e5

xmin2 = 0 * u"m"
xmax2 = 4 * u"m"
ymin2(x) = x^2 * u"m"
ymax2(x) = x^3
func2(y, x) = y * x^3

@test abs(dblquad(func2, xmin2, xmax2, ymin2, ymax2, rtol = rtol)[1] -
          res * u"m^10") < rtol * 1e2 * u"m^10"

#%%
rtol = 1e-10

xmin = 0
xmax = 4
ymin(x) = x^2
ymax(x) = x^3
zmin(x, y) = 3
zmax(x, y) = 4

func1(z, y, x) = y * x^3 * sin(z)
res=241664/5*(cos(3) - cos(4))
@test abs(tplquad(func1, xmin, xmax, ymin, ymax, zmin, zmax, rtol = rtol)[1] -res) < rtol * 1e2

@test abs(tplquad(func1, xmin, xmax, ymin, ymax, zmin, zmax, rtol = rtol, method=:vegas)[1] -res) < rtol * 1e2
@test abs(tplquad(func1, xmin, xmax, ymin, ymax, zmin, zmax, rtol = rtol, method=:divonne)[1] -res) < rtol * 1e2
@test abs(tplquad(func1, xmin, xmax, ymin, ymax, zmin, zmax, rtol = rtol, method=:suave)[1] -res) < rtol * 1e2

xmin2 = 0u"m"
xmax2 = 4u"m"
ymin2(x) = x^2*u"m"
ymax2(x) = x^3
zmin2(x, y) = 3
zmax2(x, y) = 4

func2(z, y, x) = y * x^3 * sin(z)
≈
isapprox(tplquad(func2, xmin2, xmax2, ymin2, ymax2, zmin2, zmax2, rtol = rtol)[1] , res*u"m^10", atol=rtol*1e2u"m^10")
