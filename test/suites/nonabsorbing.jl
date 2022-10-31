function test_03_bh_dielectric()
    m = 1.55
    λ0 = 0.6328
    radius = 0.525
    x = 2π * radius / λ0
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 3.10543 atol = 0.00001
    @test qsca ≈ 3.10543 atol = 0.00001
    @test qback ≈ 2.92534 atol = 0.00001
    @test g ≈ 0.63314 atol = 0.00001
end

function test_04_non_dielectric()
    m = 1.55 - im * 0.1
    λ0 = 0.6328
    radius = 0.525
    x = 2π * radius / λ0
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 2.86165188243 atol = 1e-7
    @test qsca ≈ 1.66424911991 atol = 1e-7
    @test qback ≈ 0.20599534080 atol = 1e-7
    @test g ≈ 0.80128972639 atol = 1e-7
end

function test_05_wiscombe_non_absorbing()
    # MIEV0 Test Case 5
    m = 0.75 + 0.0im
    x = 0.099
    s1 = 1.81756e-8 - 1.64810e-4 * im
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000007 atol = 1e-6
    @test g ≈ 0.001448 atol = 1e-6
    @test qback ≈ G atol = 1e-6

    # MIEV0 Test Case 6
    m = 0.75 + 0.0im
    x = 0.101
    s1 = 2.04875e-8 - im * 1.74965e-4
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000008 atol = 1e-6
    @test g ≈ 0.001507 atol = 1e-6
    @test qback ≈ G atol = 1e-6

    # MIEV0 Test Case 7
    m = 0.75 + 0.0im
    x = 10.0
    s1 = -1.07857 - im * 3.60881e-2
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.232265 atol = 1e-6
    @test g ≈ 0.896473 atol = 1e-6
    @test qback ≈ G atol = 1e-6

    # MIEV0 Test Case 8
    m = 0.75 + 0.0im
    x = 1000.0
    s1 = 1.70578e1 - im * 4.84251e2
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 1.997908 atol = 1e-6
    @test g ≈ 0.844944 atol = 1e-6
    @test qback ≈ G atol = 1e-6
end

function test_05_old_wiscombe_non_absorbing()
    # OLD MIEV0 Test Case 1
    m = 1.5 + 0.0im
    x = 10
    s1 = 4.322 - im * 4.868
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.882 atol = 1e-4
    @test qback ≈ G atol = 1e-4

    # OLD MIEV0 Test Case 2
    m = 1.5 + 0.0im
    x = 100
    s1 = 4.077e1 + im * 5.175e1
    G = abs2(2 * s1 / x)
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.0944 atol = 1e-4
    @test qback ≈ G atol = 1e-4

    # OLD MIEV0 Test Case 3
    m = 1.5 + 0.0im
    x = 1000
    G = 4 * 2.576e6 / x^2
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.0139 atol = 1e-4
    @test qback ≈ G atol = 1e-3

    # OLD MIEV0 Test Case 4
    m = 1.5 + 0.0im
    x = 5000
    G = 4 * 2.378e8 / x^2
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.0086 atol = 1e-4
    @test qback ≈ G atol = 3e-3
end

@testset "Non-Absorbing" begin
    test_03_bh_dielectric()
    test_04_non_dielectric()
    test_05_wiscombe_non_absorbing()
    test_05_old_wiscombe_non_absorbing()
end
