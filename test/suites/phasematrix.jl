function test_mie_phase_matrix_basic()
    m = 1.5 - im * 1.5
    x = 2
    μ = LinRange(-1, 1, 1000)

    p = mie_phase_matrix(m, x, μ)
    p00 = i_unpolarized(m, x, μ)

    @test all(x -> isapprox(x...), zip(p[1, 1, :], p00))
end

function test_mie_phase_matrix_mu_scalar()
    @test size(mie_phase_matrix(1.5, 2.0, 0.0)) == (4, 4)
end

function test_mie_phase_matrix_symmetry()
    p = mie_phase_matrix(1.5, 2.0, LinRange(-1, 1, 10))

    @test all(x -> isapprox(x...), zip(p[1, 2, :], p[2, 1, :]))
    @test all(x -> isapprox(x...), zip(p[3, 4, :], -p[4, 3, :]))
end

function test_mie_phase_matrix_unity()
    m = 1.5 - im * 1.5
    x = 2
    μ = LinRange(-1, 1, 1000)

    p = mie_phase_matrix(m, x, μ)

    @test all(
        x -> isapprox(x...),
        zip(p[1, 1, :] .^ 2, p[1, 2, :] .^ 2 + p[3, 3, :] .^ 2 + p[3, 4, :] .^ 2),
    )
end

function test_mie_phase_matrix_bohren()
    mm = [
        000.00 0.100000E+01 0.000000E+00 0.100000E+01 0.000000E+00
        009.00 0.785390E+00 -0.459811E-02 0.999400E+00 0.343261E-01
        018.00 0.356897E+00 -0.458541E-01 0.986022E+00 0.160184E+00
        027.00 0.766119E-01 -0.364744E+00 0.843603E+00 0.394076E+00
        036.00 0.355355E-01 -0.534997E+00 0.686967E+00 -0.491787E+00
        045.00 0.701845E-01 0.959953E-02 0.959825E+00 -0.280434E+00
        054.00 0.574313E-01 0.477927E-01 0.985371E+00 0.163584E+00
        063.00 0.219660E-01 -0.440604E+00 0.648043E+00 0.621216E+00
        072.00 0.125959E-01 -0.831996E+00 0.203255E+00 -0.516208E+00
        081.00 0.173750E-01 0.341670E-01 0.795354E+00 -0.605182E+00
        090.00 0.124601E-01 0.230462E+00 0.937497E+00 0.260742E+00
        099.00 0.679093E-02 -0.713472E+00 -0.717397E-02 0.700647E+00
        108.00 0.954239E-02 -0.756255E+00 -0.394748E-01 -0.653085E+00
        117.00 0.863419E-02 -0.281215E+00 0.536251E+00 -0.795835E+00
        126.00 0.227421E-02 -0.239612E+00 0.967602E+00 0.795798E-01
        135.00 0.543998E-02 -0.850804E+00 0.187531E+00 -0.490882E+00
        144.00 0.160243E-01 -0.706334E+00 0.495254E+00 -0.505781E+00
        153.00 0.188852E-01 -0.891081E+00 0.453277E+00 -0.226817E-01
        162.00 0.195254E-01 -0.783319E+00 -0.391613E+00 0.482752E+00
        171.00 0.301676E-01 -0.196194E+00 -0.962069E+00 0.189556E+00
        180.00 0.383189E-01 0.000000E+00 -0.100000E+01 0.000000E+00
    ]

    m = 1.55
    x = 5.213
    θ = LinRange(0, 180, 21)

    μ = cosd.(θ)
    p = mie_phase_matrix(m, x, μ, norm = :bohren)

    @test all(x -> isapprox(x...), zip(θ, mm[:, 1]))
    @test all(x -> isapprox(x...; atol=1e-3), zip(p[1, 1, :] / p[1, 1, 1], mm[:, 2]))
    @test all(x -> isapprox(x...; atol=1e-3), zip(-p[1, 2, :] ./ p[1, 1, :], mm[:, 3]))
    @test all(x -> isapprox(x...; atol=1e-3), zip(p[3, 3, :] ./ p[1, 1, :], mm[:, 4]))
    @test all(x -> isapprox(x...; atol=1e-3), zip(p[4, 3, :] ./ p[1, 1, :], mm[:, 5]))
end

@testset "MiePhaseMatrix" begin
    test_mie_phase_matrix_basic()
    test_mie_phase_matrix_mu_scalar()
    test_mie_phase_matrix_symmetry()
    test_mie_phase_matrix_unity()
    test_mie_phase_matrix_bohren()
end
