using MultiQuad, Test

@testset "quad" begin
    rtol = 1e-10

    xmin = 0
    xmax = 4
    func1(x) = x^2 * exp(-x)
    res = 2 - 26 / exp(4)
    int, err = quad(func1, xmin, xmax, rtol = rtol)
    @test isapprox(
        int,
        res,
        atol = err,
    )
end
