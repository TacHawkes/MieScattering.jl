"""
Mie scattering calculations for perfect spheres based on miepython.
Extensive documentation is at (https://miepython.readthedocs.io).

`MieScattering.jl` is a Julia package to calculate light scattering of
a plane wave by non-np.absorbing, partially-np.absorbing, or perfectly conducting
spheres.

The extinction efficiency, scattering efficiency, backscattering, and
scattering asymmetry for a sphere with complex index of refraction m,
diameter d, and wavelength lambda can be found by:

```
    qext, qsca, qback, g = ez_mie(m, d, λ0)
```
The normalized scattering values for angles µ=cos(θ) are:

```
    Ipar, Iper = ez_intensities(m, d, λ0, µ)
```

If the size parameter is known, then use:

```
    mie(m, x)
```

Mie scattering amplitudes S1 and S2 (complex numbers):

```
    mie_S1_S2(m, x, μ)
```
Normalized Mie scattering intensities for angles µ=cos(θ):

```
    i_per(m, x, µ)
    i_par(m, x, µ)
    i_unpolarized(m, x, µ)
```
"""
module MieScattering

using Polyester

include("mie.jl")
include("precompile.jl")

export  mie, mie_S1_S2, mie_cdf, mie_mu_with_uniform_cdf, generate_mie_costheta,
        i_unpolarized, i_par, i_per, ez_mie, ez_intensities, mie_phase_matrix

end
