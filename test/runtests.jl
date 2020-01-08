using MultiQuad, Test, Unitful

@testset "quad" begin
    rtol = 1e-10

    xmin = 0
    xmax = 4
    func1(x) = x^2 * exp(-x)
    @test abs(quad(func1, xmin, xmax, rtol = rtol)[1] -
              1.5237933888929114) < rtol * 1e2
    @test abs(quad(func1, xmin, xmax, rtol = rtol, method = :vegas)[1] -
              1.5237995075453157) <= rtol * 1e2
    @test abs(quad(func1, xmin, xmax, rtol = rtol, method = :suave)[1] -
              1.5288445051126394) <= rtol * 1e2

    xmin2 = 0 * u"m"
    xmax2 = 4 * u"m"
    func2(x) = x^2
    @test abs(quad(func2, xmin2, xmax2, rtol = rtol)[1] -
              21.333333333333332 * u"m^3") < rtol * 1e2 * u"m^3"
    @test abs(quad(func2, xmin2, xmax2, rtol = rtol, method = :vegas)[1] -
              21.333438438582508 * u"m^3") <= rtol * 1e2 * u"m^3"
    @test abs(quad(func2, xmin2, xmax2, rtol = rtol, method = :suave)[1] -
              21.401196340660448 * u"m^3") <= rtol * 1e2 * u"m^3"

    @test_throws ErrorException quad(func1, 0, 4, method = :none)
end

@testset "dblquad" begin
    rtol = 1e-10

    xmin = 0
    xmax = 4
    ymin(x) = x^2
    ymax(x) = x^3

    func1(y, x) = y * x^3
    @test abs(dblquad(func1, xmin, xmax, ymin, ymax, rtol = rtol)[1] -
              48332.79999999998) < rtol * 1e2
    @test abs(dblquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        rtol = rtol,
        method = :vegas,
    )[1] - 48332.56784003042) <= rtol * 1e2
    @test abs(dblquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        rtol = rtol,
        method = :suave,
    )[1] - 48317.671265074176) <= rtol * 1e2
    @test abs(dblquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        rtol = rtol,
        method = :divonne,
    )[1] - 48345.25072378812) <= rtol * 1e2

    xmin2 = 0 * u"m"
    xmax2 = 4 * u"m"
    ymin2(x) = x^2 * u"m"
    ymax2(x) = x^3
    func2(y, x) = y * x^3

    @test abs(dblquad(func2, xmin2, xmax2, ymin2, ymax2, rtol = rtol)[1] -
              48332.79999999998 * u"m^10") < rtol * 1e2 * u"m^10"
    @test abs(dblquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        rtol = rtol,
        method = :vegas,
    )[1] - 48332.56784003042 * u"m^10") <= rtol * 1e2 * u"m^10"
    @test abs(dblquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        rtol = rtol,
        method = :suave,
    )[1] - 48317.671265074176 * u"m^10") <= rtol * 1e2 * u"m^10"
    @test abs(dblquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        rtol = rtol,
        method = :divonne,
    )[1] - 48345.25072378812 * u"m^10") <= rtol * 1e2 * u"m^10"

    @test_throws ErrorException dblquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        rtol = rtol,
        method = :none,
    )
end

@testset "tblquad" begin
    rtol = 1e-10

    xmin = 0
    xmax = 4
    ymin(x) = x^2
    ymax(x) = x^3
    zmin(x, y) = 3
    zmax(x, y) = 4

    func1(z, y, x) = y * x^3 * sin(z)
    @test abs(tplquad(func1, xmin, xmax, ymin, ymax, zmin, zmax, rtol = rtol)[1] + 16256.682941203506) < rtol * 1e2
    @test abs(tplquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
        rtol = rtol,
        method = :vegas,
    )[1] + 16256.948833335844) < rtol * 1e2
    @test abs(tplquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
        rtol = rtol,
        method = :suave,
    )[1] + 16254.99067642033) < rtol * 1e2
    @test abs(tplquad(
        func1,
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
        rtol = rtol,
        method = :divonne,
    )[1] + 16262.818635298641) < rtol * 1e2

    xmin2 = 0u"m"
    xmax2 = 4u"m"
    ymin2(x) = x^2*u"m"
    ymax2(x) = x^3
    zmin2(x, y) = 3
    zmax2(x, y) = 4

    func2(z, y, x) = y * x^3 * sin(z)
    @test abs(tplquad(func2, xmin2, xmax2, ymin2, ymax2, zmin2, zmax2, rtol = rtol)[1] + 16256.682941203506u"m^10") < rtol * 1e2u"m^10"
    @test abs(tplquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        zmin2,
        zmax2,
        rtol = rtol,
        method = :vegas,
    )[1] + 16256.948833335844u"m^10") < rtol * 1e2u"m^10"
    @test abs(tplquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        zmin2,
        zmax2,
        rtol = rtol,
        method = :suave,
    )[1] + 16254.99067642033u"m^10") < rtol * 1e2u"m^10"
    @test abs(tplquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        zmin2,
        zmax2,
        rtol = rtol,
        method = :divonne,
    )[1] + 16262.818635298641u"m^10") < rtol * 1e2u"m^10"

    @test_throws ErrorException tplquad(
        func2,
        xmin2,
        xmax2,
        ymin2,
        ymax2,
        zmin2,
        zmin2,
        rtol = rtol,
        method = :none,
    )
end
