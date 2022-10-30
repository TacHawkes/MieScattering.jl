using MieScattering, GLMakie

x = LinRange(0.1, 100, 300)
qext, qsca, qback, g = mie(1.5, x)

fig = Figure()
ax = fig[1,1] = Axis(fig, xlabel="Size parameter (-)", ylabel="Qext", title="Comparison of extinction for absorbing and non-absorbing spheres")
lines!(x, qext, color=:red, label="1.5")

qext, qsca, qback, g = mie(1.5 -0.1*im, x)
lines!(x, qext, color=:blue, label="1.5-0.1j")

axislegend()
fig
