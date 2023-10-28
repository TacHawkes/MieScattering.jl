function test_12_scatter_function()
    x = 1.0
    m = 1.5 - im * 1.0
    θ = 0:30:180
    μ = cosd.(θ)

    qext, qsca, qback, g = mie(m, x)
    S1, S2 = mie_S1_S2(m, x, μ)
    S1 *= √(π * x^2 * qext)
    S2 *= √(π * x^2 * qext)

    @test real(S1[1]) ≈ 0.584080 atol = 1e-6
    @test imag(S1[1]) ≈ 0.190515 atol = 1e-6
    @test real(S2[1]) ≈ 0.584080 atol = 1e-6
    @test imag(S2[1]) ≈ 0.190515 atol = 1e-6

    @test real(S1[2]) ≈ 0.565702 atol = 1e-6
    @test imag(S1[2]) ≈ 0.187200 atol = 1e-6
    @test real(S2[2]) ≈ 0.500161 atol = 1e-6
    @test imag(S2[2]) ≈ 0.145611 atol = 1e-6

    @test real(S1[3]) ≈ 0.517525 atol = 1e-6
    @test imag(S1[3]) ≈ 0.178443 atol = 1e-6
    @test real(S2[3]) ≈ 0.287964 atol = 1e-6
    @test imag(S2[3]) ≈ 0.041054 atol = 1e-6

    @test real(S1[4]) ≈ 0.456340 atol = 1e-6
    @test imag(S1[4]) ≈ 0.167167 atol = 1e-6
    @test real(S2[4]) ≈ 0.0362285 atol = 1e-6
    @test imag(S2[4]) ≈ -0.0618265 atol = 1e-6

    @test real(S1[5]) ≈ 0.400212 atol = 1e-6
    @test imag(S1[5]) ≈ 0.156643 atol = 1e-6
    @test real(S2[5]) ≈ -0.174875 atol = 1e-6
    @test imag(S2[5]) ≈ -0.122959 atol = 1e-6

    @test real(S1[6]) ≈ 0.362157 atol = 1e-6
    @test imag(S1[6]) ≈ 0.149391 atol = 1e-6
    @test real(S2[6]) ≈ -0.305682 atol = 1e-6
    @test imag(S2[6]) ≈ -0.143846 atol = 1e-6

    @test real(S1[7]) ≈ 0.348844 atol = 1e-6
    @test imag(S1[7]) ≈ 0.146829 atol = 1e-6
    @test real(S2[7]) ≈ -0.348844 atol = 1e-6
    @test imag(S2[7]) ≈ -0.146829 atol = 1e-6
end

function test_13_unity_normalization()
    x = 1.0
    m = 1.5 - im * 1.0
    θ = 0:30:180
    μ = cosd.(θ)
    qext, qsca, qback, g = mie(m, x)
    S1, S2 = mie_S1_S2(m, x, μ, norm = :wiscombe)

    @test real(S1[1]) ≈ 0.584080 atol = 1e-6
    @test imag(S1[1]) ≈ 0.190515 atol = 1e-6
    @test real(S2[1]) ≈ 0.584080 atol = 1e-6
    @test imag(S2[1]) ≈ 0.190515 atol = 1e-6

    @test real(S1[2]) ≈ 0.565702 atol = 1e-6
    @test imag(S1[2]) ≈ 0.187200 atol = 1e-6
    @test real(S2[2]) ≈ 0.500161 atol = 1e-6
    @test imag(S2[2]) ≈ 0.145611 atol = 1e-6

    @test real(S1[3]) ≈ 0.517525 atol = 1e-6
    @test imag(S1[3]) ≈ 0.178443 atol = 1e-6
    @test real(S2[3]) ≈ 0.287964 atol = 1e-6
    @test imag(S2[3]) ≈ 0.041054 atol = 1e-6

    @test real(S1[4]) ≈ 0.456340 atol = 1e-6
    @test imag(S1[4]) ≈ 0.167167 atol = 1e-6
    @test real(S2[4]) ≈ 0.0362285 atol = 1e-6
    @test imag(S2[4]) ≈ -0.0618265 atol = 1e-6

    @test real(S1[5]) ≈ 0.400212 atol = 1e-6
    @test imag(S1[5]) ≈ 0.156643 atol = 1e-6
    @test real(S2[5]) ≈ -0.174875 atol = 1e-6
    @test imag(S2[5]) ≈ -0.122959 atol = 1e-6

    @test real(S1[6]) ≈ 0.362157 atol = 1e-6
    @test imag(S1[6]) ≈ 0.149391 atol = 1e-6
    @test real(S2[6]) ≈ -0.305682 atol = 1e-6
    @test imag(S2[6]) ≈ -0.143846 atol = 1e-6

    @test real(S1[7]) ≈ 0.348844 atol = 1e-6
    @test imag(S1[7]) ≈ 0.146829 atol = 1e-6
    @test real(S2[7]) ≈ -0.348844 atol = 1e-6
    @test imag(S2[7]) ≈ -0.146829 atol = 1e-6
end

function test_i_unpolarized_01()
    m = 1.5 - im * 1.5
    x = 2
    μ = LinRange(-1, 1, 1000)
    qext, qsca, _, _ = mie(m, x)
    expected = [qsca / qext, 1.0, 4π, qsca, qext, qsca * 4 * π * x^2, qsca * π * x^2]

    for (i, norm) in enumerate([:albedo, :one, :four_pi, :qsca, :qext, :bohren, :wiscombe])
        intensity = i_unpolarized(m, x, μ; norm)
        total = 2π * (μ[2] - μ[1]) * sum(intensity)
        @test total / expected[i] ≈ 1.0 atol = 4e-3
    end
end

function test_i_par_i_per_01()
    m = 1.5 - im * 1.5
    x = 2
    μ = LinRange(-1, 1, 1000)
    qext, qsca, _, _ = mie(m, x)
    expected = [qsca / qext, 1.0, 4π, qsca, qext, qsca * 4 * π * x^2, qsca * π * x^2]
    for (i, norm) in enumerate([:albedo, :one, :four_pi, :qsca, :qext, :bohren, :wiscombe])
        iper = i_per(m, x, μ; norm)
        total1 = 2π * (μ[2] - μ[1]) * sum(iper)
        ipar = i_par(m, x, μ; norm)
        total2 = 2π * (μ[2] - μ[1]) * sum(ipar)
        total = (total1 + total2) / 2
        @test total / expected[i] ≈ 1.0 atol = 5e-3
    end
end

function test_molecular_hydrogen()
    m = 1.00013626
    x = 0.0006403246172921872
    mu = LinRange(-1, 1, 100)
    ph = i_unpolarized(m, x, mu)
    @test ph[2] ≈ 0.1169791 atol=1e-5
end

@testset "AngleScattering" begin
    test_12_scatter_function()
    test_13_unity_normalization()
    test_i_unpolarized_01()
    test_i_par_i_per_01()
    test_molecular_hydrogen()
end
