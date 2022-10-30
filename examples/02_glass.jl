using MieScattering, GLMakie

radius = 2 # in microns
λ0 = LinRange(0.2, 1.2, 200)
x = 2π*radius./λ0

# from https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
m2 = @. 1 + 1.03961212 / (1 - 0.00600069867 / λ0^2)
m2 += @. 0.231792344 / (1 - 0.0200179144 / λ0^2)
m2 += @. 1.01046945 / (1 - 103.560653 / λ0^2)
m = .√(m2)
qext, qsca, qback, g = mie(m, x)

fig = Figure()
ax = fig[1,1] = Axis(fig, xlabel="Wavelength [nm]", ylabel="Scattering efficiency", title="BK7 glass spheres 4 micron diameter")
lines!(λ0*1000, qsca)

fig
