function test_06_wiscombe_water_absorbing()
    # MIEV0 Test Case 9
    m = 1.33 - im * 0.00001
    x = 1.0
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.093923 atol = 1e-6
    @test g ≈ 0.184517 atol = 1e-6

    # MIEV0 Test Case 10
    m = 1.33 - im * 0.00001
    x = 100.0
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.096594 atol = 1e-6
    @test g ≈ 0.868959 atol = 1e-6

    # MIEV0 Test Case 11
    m = 1.33 - im * 0.00001
    x = 10000.0
    qext, qsca, qback, g = mie(m, x)
    @test g ≈ 0.907840 atol = 1e-6
    @test qsca ≈ 1.723857 atol = 1e-6
end

function test_07_wiscombe_absorbing()
    # MIEV0 Test Case 12
    m = 1.5 - im
    x = 0.055
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000011 atol = 1e-6
    @test g ≈ 0.000491 atol = 1e-6

    # MIEV0 Test Case 13
    m = 1.5 - im
    x = 0.056
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.000012 atol = 1e-6
    @test g ≈ 0.000509 atol = 1e-6

    # MIEV0 Test Case 14
    m = 1.5 - im
    x = 1
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 0.6634538 atol = 1e-6
    @test g ≈ 0.192136 atol = 1e-6

    # MIEV0 Test Case 15
    m = 1.5 - im
    x = 100.0
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 1.283697 atol = 1e-3
    @test qext ≈ 2.097502 atol = 1e-2
    @test g ≈ 0.850252 atol = 1e-3

    # MIEV0 Test Case 16
    m = 1.5 - im
    x = 10000
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 1.236575 atol = 1e-6
    @test qext ≈ 2.004368 atol = 1e-6
    @test g ≈ 0.846309 atol = 1e-6
end

function test_08_wiscombe_more_absorbing()
    # MIEV0 Test Case 17
    m = 10.0 - im * 10.0
    x = 1
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 2.049405 atol = 1e-6
    @test g ≈ -0.110664 atol = 1e-6

    # MIEV0 Test Case 18
    m = 10.0 - im * 10.0
    x = 100
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 1.836785 atol = 1e-6
    @test g ≈ 0.556215 atol = 1e-6

    # MIEV0 Test Case 19
    m = 10.0 - im * 10.0
    x = 10000
    qext, qsca, qback, g = mie(m, x)
    @test qsca ≈ 1.795393 atol = 1e-6
    @test g ≈ 0.548194 atol = 1e-6
end

function test_09_single_nonmagnetic()
    m = 1.5 - im * 0.5
    x = 2.5
    qext, qsca, qback, g = mie(m, x)
    @test qext ≈ 2.562873497454734 atol = 1e-7
    @test qsca ≈ 1.097071819088392 atol = 1e-7
    @test qback ≈ 0.123586468179818 atol = 1e-7
    @test g ≈ 0.748905978948507 atol = 1e-7
end

@testset "Absorbing" begin
    test_06_wiscombe_water_absorbing()
    test_07_wiscombe_absorbing()
    test_08_wiscombe_more_absorbing()
    test_09_single_nonmagnetic()
end
