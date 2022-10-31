"""
Plot the scattering cross section for 100nm gold spheres.
The resulting graph is as a function of wavelength.
"""

using MieScattering
using CairoMakie

# from https://refractiveindex.info/?shelf = main&book = Au&page = Johnson
# wavelength in microns
ref_lam = [ 0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164,
            0.2214, 0.2262, 0.2313, 0.2371, 0.2426, 0.2490, 0.2551, 0.2616,
            0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107, 0.3204, 0.3315,
            0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509,
            0.4714, 0.4959, 0.5209, 0.5486, 0.5821, 0.6168, 0.6595, 0.7045,
            0.7560, 0.8211, 0.8920, 0.9840, 1.0880, 1.2160, 1.3930, 1.6100, 1.9370]

ref_n = [   1.28, 1.32, 1.34, 1.33, 1.33, 1.30, 1.30, 1.30, 1.30, 1.31, 1.30,
            1.32, 1.32, 1.33, 1.33, 1.35, 1.38, 1.43, 1.47, 1.49, 1.53, 1.53,
            1.54, 1.48, 1.48, 1.50, 1.48, 1.46, 1.47, 1.46, 1.45, 1.38, 1.31,
            1.04, 0.62, 0.43, 0.29, 0.21, 0.14, 0.13, 0.14, 0.16, 0.17, 0.22,
            0.27, 0.35, 0.43, 0.56, 0.92]

ref_k = [   1.188, 1.203, 1.226, 1.251, 1.277, 1.304, 1.350, 1.387, 1.427,
            1.460, 1.497, 1.536, 1.577, 1.631, 1.688, 1.749, 1.803, 1.847,
            1.869, 1.878, 1.889, 1.893, 1.898, 1.883, 1.871, 1.866, 1.895,
            1.933, 1.952, 1.958, 1.948, 1.914, 1.849, 1.833, 2.081, 2.455,
            2.863, 3.272, 3.697, 4.103, 4.542, 5.083, 5.663, 6.350, 7.150,
            8.145, 9.519, 11.21, 13.78]

radius = 0.1 # in microns
m = ref_n - im * ref_k
x = 2π * radius ./ ref_lam
cross_section_area = π * radius ^ 2
mu_a = 4 * π * ref_k ./ ref_lam    # nm
qext, qsca, qback, g = mie(m, x)

sca_cross_section = qsca * cross_section_area
abs_cross_section = (qext - qsca) * cross_section_area

fig = Figure()

ax = fig[1,1] = Axis(fig, title="Gold Spheres 200nm diameter")
hidexdecorations!(ax)
scatter!(ref_lam * 1000, ref_n, color=:blue, marker=:circle)
scatter!(ref_lam * 1000, -ref_k, color=:red, marker=:rect)
text!(700, 1; text="real refractive index", color=:blue)
text!(1100, -6; text="imaginary refractive index", color=:red)
ylims!(-14,5)

ax = fig[2,1] = Axis(fig, ylabel="Absorption Depth [nm]")
scatter!(ref_lam * 1000, 1000 ./ mu_a, color=:blue, marker=:circle)

ax = fig[3,1] = Axis(fig, xlabel="Wavelength [nm]", ylabel="Cross Section [µm²]")
scatter!(ref_lam * 1000, abs_cross_section, color=:blue, marker=:circle)
scatter!(ref_lam * 1000, sca_cross_section, color=:red, marker=:rect)
text!(700, 0.01; text="absorption", color=:blue)
text!(750, 0.1; text="scattering", color=:red)

fig