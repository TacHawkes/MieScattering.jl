using MieScattering, GLMakie, BenchmarkTools

num = 100
radius = 0.5 # in microns
λ0 = LinRange(0.2, 1.2, 200)
x = 2π * radius ./ λ0

# from https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon - 24.0C
m2 = @. 1.0 + 5.666959820E-1 / (1.0 - 5.084151894E-3 / λ0^2)
m2 += @. 1.731900098E-1 / (1.0 - 1.818488474E-2 / λ0^2)
m2 += @. 2.095951857E-2 / (1.0 - 2.625439472E-2 / λ0^2)
m2 += @. 1.125228406E-1 / (1.0 - 1.073842352E1 / λ0^2)
m = .√(m2)
qext, qsca, qback, g = mie(m, x)

fig = Figure()
ax =
    fig[1, 1] = Axis(
        fig,
        xlabel = "Wavelength [nm]",
        ylabel = "Scattering Cross Section (µm^2)",
        title = "Water Droplets (1 µm diameter)",
    )
lines!(λ0 * 1000, qsca)

fig
