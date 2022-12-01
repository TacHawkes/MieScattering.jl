using SnoopPrecompile

@precompile_setup begin
    m_sphere = 1.0
    n_water = 4 / 3
    d = 1000
    λ0 = LinRange(300, 800, 50)

    x = 1.0
    m = 1.5 - im * 1.0
    θ = 0:30:180
    μ = cosd.(θ)
    @precompile_all_calls begin
        qext, qsca, qback, g = ez_mie(m_sphere, d, λ0, n_water)

        qext, qsca, qback, g = mie(m, x)
        S1, S2 = mie_S1_S2(m, x, μ)
    end
end
