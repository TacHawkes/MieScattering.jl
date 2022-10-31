function test_10_small_spheres()
    # MIEV0 Test Case 5
    m = 0.75
    x = 0.099

    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000007 atol = 1e-6
    @test g ≈ 0.001448 atol = 1e-6

    # MIEV0 Test Case 6
    m = 0.75
    x = 0.101

    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000008 atol = 1e-6
    @test g ≈ 0.001507 atol = 1e-6

    m = 1.5 - im
    x = 0.055
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.101491 atol = 1e-6
    @test g ≈ 0.000491 atol = 1e-6
    x = 0.056
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.103347 atol = 1e-6
    @test g ≈ 0.000509 atol = 1e-6

    m = 1e-10 - im * 1e10
    x = 0.099
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000321 atol = 1e-6
    @test g ≈ -0.397357 atol = 1e-4
    x = 0.101
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000348 atol = 1e-6
    @test g ≈ -0.397262 atol = 1e-6

    m = 0 - im * 1e10
    x = 0.099
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000321 atol = 1e-6
    @test g ≈ -0.397357 atol = 1e-4
    x = 0.101
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 0.000348 atol = 1e-6
    @test g ≈ -0.397262 atol = 1e-4
end

@testset "Small" begin
    test_10_small_spheres()
end
