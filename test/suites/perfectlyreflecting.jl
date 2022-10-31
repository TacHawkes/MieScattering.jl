function test_11_wiscombe_perfectly_reflecting()
    # MIEV0 Test Case 0
    m = 0
    x = 0.001

    qext, qsca, qback, g = mie(m, x)
    @test isapprox(qsca, 3.3333e-12, atol = 1e-13)

    # MIEV0 Test Case 1
    m = 0
    x = 0.099

    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000321 atol = 1e-4
    @test g ≈ -0.397357 atol = 1e-3

    # MIEV0 Test Case 2
    m = 0
    x = 0.101

    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000348 atol = 1e-6
    @test g ≈ -0.397262 atol = 1e-6

    # MIEV0 Test Case 3
    m = 0
    x = 100

    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.008102 atol = 1e-6
    @test g ≈ 0.500926 atol = 1e-6

    # MIEV0 Test Case 4
    m = 0
    x = 10000

    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.000289 atol = 1e-6
    @test g ≈ 0.500070 atol = 1e-6
end

@testset "Perfectly reflecting" begin
    test_11_wiscombe_perfectly_reflecting()
end
