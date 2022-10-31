function test_nb1_x()
    N = 500
    m = 1.5
    x = LinRange(0.1, 20, N)
    qext, qsca, qback, g = mie(m, x)

    @test qsca[1] ≈ 2.3084093592198083e-05 atol=1e-6
    @test qsca[101] ≈ 4.105960809066763 atol=1e-6
    @test qsca[201] ≈ 1.9947867190110644 atol=1e-6
    @test qsca[301] ≈ 2.4652591512196405 atol=1e-6
    @test qsca[401] ≈ 2.472171798724846 atol=1e-6
    @test qsca[500] ≈ 2.03583698038088 atol=1e-6
end

function test_nb1_rho()
    N = 500
    m = 1.5
    ρ = LinRange(0.1, 20, N)
    x15 = ρ/2/(m-1)
    qext, scal5, qback, g = mie(m, x15)

    m = 1.1
    x11 = ρ/2/(m-1)
    qext, scal1, qback, g = mie(m, x11)

    @test scal1[1] ≈ 0.0006616369953521216 atol=1e-6
    @test scal1[100] ≈ 3.449616595439377 atol=1e-6
    @test scal1[200] ≈ 1.6837703285684387 atol=1e-6
    @test scal1[300] ≈ 2.3167184401740495 atol=1e-6
    @test scal1[400] ≈ 2.218210809017406 atol=1e-6
    @test scal1[500] ≈ 1.876467571615533 atol=1e-6

    @test scal5[1] ≈ 2.3084093592198083e-05 atol=1e-6
    @test scal5[100] ≈ 4.07295075914037 atol=1e-6
    @test scal5[200] ≈ 1.8857586341949146 atol=1e-6
    @test scal5[300] ≈ 2.464763930426085 atol=1e-6
    @test scal5[400] ≈ 2.430569030744473 atol=1e-6
    @test scal5[500] ≈ 2.03583698038088 atol=1e-6
end

function test_nb1_spheres()
    N = 500
    m = 1.0
    r = 500 # nm
    λ0 = LinRange(300, 800, N)

    mwater = 4/3
    mm = m / mwater
    xx = 2π*r*mwater ./ λ0

    qext, qsca, qback, g = mie(mm, xx)

    @test qsca[1] ≈ 1.5525047718022498 atol=1e-6
    @test qsca[100] ≈ 2.1459528526672678 atol=1e-6
    @test qsca[200] ≈ 2.365171370327149 atol=1e-6
    @test qsca[300] ≈ 2.2039860928542128 atol=1e-6
    @test qsca[400] ≈ 1.9261758931397088 atol=1e-6
    @test qsca[500] ≈ 1.640006561518987 atol=1e-6
end

function test_nb1_ezmie()
    m_sphere = 1.0
    n_water = 4/3
    d = 1000
    λ0 = LinRange(300, 800, 50)
    qext, qsca, qback, g = ez_mie(m_sphere, d, λ0, n_water)

    @test qsca[1] ≈ 1.5525047718022498 atol=1e-6
    @test qsca[10] ≈ 2.107970892634116 atol=1e-6
    @test qsca[20] ≈ 2.3654333205160074 atol=1e-6
    @test qsca[30] ≈ 2.213262310704816 atol=1e-6
    @test qsca[40] ≈ 1.9314911518355427 atol=1e-6
    @test qsca[50] ≈ 1.640006561518987 atol=1e-6
end

@testset "Notebook" begin
    test_nb1_x()
    test_nb1_rho()
    test_nb1_spheres()
    test_nb1_ezmie()
end
