"""
    mie(m, x)

Calculate the efficiencies for a sphere where m or x may be vectors.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere

# Output
- `qext`: the total extinction efficiency
- `qsca`: the scattering efficiency
- `qback`: the backscatter efficiency
- `g`: the average cosine of the scattering phase function
"""
function mie end

function mie(m::T, x::V; use_threads = true) where {T<:Number,V<:Real}
    TT = float(real(promote_type(T, V)))
    if iszero(real(m)) && x < TT(0.1)
        qext, qsca, qback, g = small_conducting_mie(x)
    elseif !iszero(real(m)) && abs(m) * x < TT(0.1)
        qext, qsca, qback, g = small_mie(m, x)
    else
        a, b = mie_An_Bn(m, x)

        nmax = length(a)
        n = 1:nmax

        qext = 2 * sum(i -> (2i + 1) * (real(a[i]) + real(b[i])), n) / x^2
        qsca = qext

        if !iszero(imag(m))
            qsca = 2 * sum(i -> (2i + 1) * (abs2(a[i]) + abs2(b[i])), n) / x^2
        end

        qback = abs2(sum(i -> (-1)^i * ((2i + 1)) * (a[i] - b[i]), n)) / x^2

        g = zero(TT)
        for i = 1:(nmax-1)
            asy1 = (i * (i + 2) / (i + 1)) * real(a[i] * conj(a[i+1]) + b[i] * conj(b[i+1]))
            asy2 = (2 * i + 1) / i / (i + 1) * real(a[i] * conj(b[i]))
            g += 4 * (asy1 + asy2) / qsca / x^2
        end
    end

    return qext, qsca, qback, g
end

function mie(
    m::AbstractVector{T},
    x::AbstractVector{V};
    use_threads = true,
) where {T<:Number,V<:Real}
    mlen, xlen = length(m), length(x)

    if xlen > 1 && mlen > 1 && xlen != mlen
        error("Vectors m and x passed to mie function must have the same length")
    end
    O = float(real(promote_type(T, V)))
    len = max(xlen, mlen)
    # allocate outputs
    qext = Vector{O}(undef, len)
    qsca = Vector{O}(undef, len)
    qback = Vector{O}(undef, len)
    g = Vector{O}(undef, len)

    @batch minbatch = (use_threads ? 1 : typemax(Int)) for i in eachindex(g)
        mm = mlen > 1 ? m[i] : first(m)
        xx = xlen > 1 ? x[i] : first(x)
        @inbounds qext[i], qsca[i], qback[i], g[i] = mie(mm, xx)
    end

    return qext, qsca, qback, g
end
mie(m::Number, x::AbstractVector; use_threads = true) = mie([m], x; use_threads)
mie(m::AbstractVector, x::Number; use_threads = true) = mie(m, [x]; use_threads)

"""
    small_mie(m, x)

Calculate the efficiencies for a small sphere.
Typically used for small spheres where `x<0.1`

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere

# Output
- `qext`: the total extinction efficiency
- `qsca`: the scattering efficiency
- `qback`: the backscatter efficiency
- `g`: the average cosine of the scattering phase function
"""
function small_mie(m, x)
    m2 = m^2
    x2 = x^2

    D = complex(m2 + 2 + (1 - 0.7 * m2) * x2)
    D -= (8 * m^4 - 385 * m2 + 350) * x^4 / 1400.0
    D += 2 * im * (m2 - 1) * x^3 * (1 - 0.1 * x2) / 3
    ahat1 = 2im * (m2 - 1) / 3 * (1 - 0.1 * x2 + (4 * m2 + 5) * x^4 / 1400) / D
    bhat1 = 1im * x2 * (m2 - 1) / 45 * (1 + (2 * m2 - 5) / 70 * x2)
    bhat1 /= 1 - (2 * m2 - 5) / 30 * x2
    ahat2 = im * x2 * (m2 - 1) / 15 * (1 - x2 / 14)
    ahat2 /= 2 * m2 + 3 - (2 * m2 - 7) / 14 * x2

    T = abs2(ahat1) + abs2(bhat1) + 5 / 3 * abs2(ahat2)
    temp = ahat2 + bhat1
    g = real(ahat1 * conj(temp)) / T

    qsca = 6 * x^4 * T

    if iszero(real(m))
        qext = qsca
    else
        qext = 6 * x * real(ahat1 + bhat1 + 5 * ahat2 / 3)
    end

    sback = 1.5 * x^3 * (ahat1 - bhat1 - 5 * ahat2 / 3)
    qback = 4 * abs2(sback) / x2

    return qext, qsca, qback, g
end

"""
    small_conducting_mie(m, x)

Calculate the efficiencies for a small conducting spheres.
Typically used for small conducting spheres where `x < 0.1` and
`real(m) == 0`.

# Parameters
- `x`: the size parameter of the sphere

# Output
- `qext`: the total extinction efficiency
- `qsca`: the scattering efficiency
- `qback`: the backscatter efficiency
- `g`: the average cosine of the scattering phase function
"""
function small_conducting_mie(x)
    ahat1 = (im * 2 / 3 * (1 - 0.2 * x^2)) / (1 - 0.5 * x^2 + im * 2 / 3 * x^3)
    bhat1 = (im * (x^2 - 10.0) / 30.0) / (1 + 0.5 * x^2 - im * x^3 / 3.0)
    ahat2 = im * x^2 / 30
    bhat2 = -im * x^2 / 45

    qsca = x^4 * (6 * abs2(ahat1) + 6 * abs2(bhat1) + 10 * abs2(ahat2) + 10 * abs2(bhat2))
    qext = qsca
    g = imag(ahat1) * (imag(ahat2) + imag(bhat1))
    g += imag(bhat2) * (5 / 9 * imag(ahat2) + imag(bhat1))
    g += real(ahat1) * real(bhat1)
    g *= 6 * x^4 / qsca

    qback = 9 * x^4 * abs(ahat1 - bhat1 - 5 / 3 * (ahat2 - bhat2))^2

    return qext, qsca, qback, g
end

"""
    mie_An_Bn(m, x)

Compute arrays of Mie coefficients A and B for a sphere.
This estimates the size of the arrays based on Wiscombe's formula. The length
of the arrays is chosen so that the error when the series are summed is around 1e-6.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere

# Output
- `An`, `Bn`: arrays of Mie coefficents
"""
function mie_An_Bn(m, x)
    nstop = floor(Int, x + 4.05 * x^0.33333 + 2.0) + 1
    a = Vector{ComplexF64}(undef, nstop - 1)
    b = Vector{ComplexF64}(undef, nstop - 1)

    psi_nm1 = sin(x)
    psi_n = psi_nm1 / x - cos(x)
    xi_nm1 = psi_nm1 + im * cos(x)
    xi_n = psi_n + im * (cos(x) / x + sin(x))

    if !iszero(real(m))
        D = D_calc(m, x, nstop + 1)

        for n = 1:(nstop-1)
            temp = D[n] / m + n / x
            a[n] = (temp * psi_n - psi_nm1) / (temp * xi_n - xi_nm1)
            temp = D[n] * m + n / x
            b[n] = (temp * psi_n - psi_nm1) / (temp * xi_n - xi_nm1)
            xi = (2 * n + 1) * xi_n / x - xi_nm1
            xi_nm1 = xi_n
            xi_n = xi
            psi_nm1 = psi_n
            psi_n = real(xi_n)
        end
    else
        for n = 1:(nstop-1)
            a[n] = (n * psi_n / x - psi_nm1) / (n * xi_n / x - xi_nm1)
            b[n] = psi_n / xi_n
            xi = (2 * n + 1) * xi_n / x - xi_nm1
            xi_nm1 = xi_n
            xi_n = xi
            psi_nm1 = psi_n
            psi_n = real(xi_n)
        end
    end

    return a, b
end

"""
    D_calc(m, x, N)

Compute the logarithmic derivative using best method.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `N`: order of Ricatti-Bessel function

# Output
The values of the Ricatti-Bessel function for orders from 0 to N.
"""
function D_calc(m, x, N)
    n = real(m)
    κ = abs(imag(m))
    D = Vector{ComplexF64}(undef, N)

    if n < 1 || n > 10 || κ > 10 || (x * κ) >= (3.9 - 10.8 * n + 13.78 * n^2)
        D_downwards!(m * x, N, D)
    else
        D_upwards!(m * x, N, D)
    end

    return D
end

"""
    D_upwards!(z, N, D)

Compute the logarithmic derivative by upwards recurrence.

# Parameters
- `z`: function argument
- `N`: order of Ricatti-Bessel function
- `D`: gets filled with the Ricatti-Bessel function values for orders
       from 0 to N for an argument z using the upwards recurrence relations.
"""
function D_upwards!(z, N, D)
    _exp = exp(-2im * z)
    @inbounds D[1] = -1 / z + (1 - _exp) / ((1 - _exp) / z - im * (1 + _exp))
    @inbounds for n = 2:N
        D[n] = 1 / (n / z - D[n-1]) - n / z
    end
end

"""
    D_downwards!(z, N, D)

Compute the logarithmic derivative by downwards recurrence.

# Parameters
- `z`: function argument
- `N`: order of Ricatti-Bessel function
- `D`: gets filled with the Ricatti-Bessel function values for orders
       from 0 to N for an argument z using the downwards recurrence relations.
"""
function D_downwards!(z, N, D)
    last_D = Lentz_Dn(z, N)
    @inbounds for n = N:-1:2
        last_D = n / z - 1.0 / (last_D + n / z)
        D[n-1] = last_D
    end
end

"""
    Lentz_Dn(z, N)

Compute the logarithmic derivative of the Ricatti-Bessel function.

# Parameters
- `z`: function argument
- `N`: order of Ricatti-Bessel function


# Output

This returns the Ricatti-Bessel function of order N with argument z
using the continued fraction technique of Lentz, Appl. Opt., 15,
668-671, (1976).
"""
function Lentz_Dn(z, N)
    zinv = 2.0 / z
    α = (N + 0.5) * zinv
    aj = -(N + 1.5) * zinv
    α_j1 = aj + 1 / α
    α_j2 = aj
    ratio = α_j1 / α_j2
    runratio = α * ratio

    while abs(abs(ratio) - 1.0) > 1e-12
        aj = zinv - aj
        α_j1 = 1.0 / α_j1 + aj
        α_j2 = 1.0 / α_j2 + aj
        ratio = α_j1 / α_j2
        zinv *= -1
        runratio = ratio * runratio
    end

    return -N / z + runratio
end

"""
    small_mie_conducting_S1_S2(m, x, μ)

Calculate the scattering amplitudes for small conducting spheres.
The spheres are small perfectly conducting (reflecting) spheres (`x<0.1`).
The amplitude functions have been normalized so that when integrated
over all `4𝜋` solid angles, the integral will be qext(`𝜋x²`).
The units are weird, ``sr^{-0.5}``.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(``θ``), to calculate scattering amplitudes

# Output
`S1`, `S2`: the scattering amplitudes at each angle µ [``sr^{-0.5}``]
"""
function small_mie_conducting_S1_S2(m, x, μ)
    ahat1 = 2im / 3 * (1 - 0.2 * x^2) / (1 - 0.5 * x^2 + 2im / 3 * x^3)
    bhat1 = 1im / 3 * (0.1 * x^2 - 1) / (1 + 0.5 * x^2 - 1im / 3 * x^3)
    ahat2 = 1im / 30 * x^2
    bhat2 = -1im * x^2 / 45

    S1 = @. 1.5 *
       x^3 *
       (ahat1 + bhat1 * μ + 5 / 3 * ahat2 * μ + 5 / 3 * bhat2 * (2 * μ^2 - 1))
    S2 = @. 1.5 *
       x^3 *
       (bhat1 + ahat1 * μ + 5 / 3 * bhat2 * μ + 5 / 3 * ahat2 * (2 * μ^2 - 1))

    qext = x^4 * (6 * abs2(ahat1) + 6 * abs(bhat1) + 10 * abs2(ahat2) + 10 * abs2(bhat2))

    norm = √(qext * π * x^2)
    S1 ./= norm
    S2 ./= norm

    return S1, S2
end

"""
Calculate the scattering amplitude functions for small spheres (`x<0.1`).
The amplitude functions have been normalized so that when integrated
over all `4*π` solid angles, the integral will be `qext*pi*x^2`.
The units are weird, ``sr^{-0.5}``

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(``θ``), to calculate scattering amplitudes

# Output
`S1`, `S2`: the scattering amplitudes at each angle µ [``sr^{-0.5}``]
"""
function small_mie_S1_S2(m, x, μ)
    m2 = m^2
    m4 = m2^2
    x2 = x^2
    x3 = x2 * x
    x4 = x2^2

    D = m2 + 2 + (1 - 0.7 * m2) * x2
    D -= (8 * m4 - 385 * m2 + 350) * x4 / 1400
    D += 2im * (m2 - 1) * x3 * (1 - 0.1 * x2) / 3

    ahat1 = 2im * (m2 - 1) / 3 * (1 - 0.1 * x2 + (4 * m2 + 5) * x4 / 1400) / D

    bhat1 = im * x2 * (m2 - 1) / 45 * (1 + (2 * m2 - 5) / 70 * x2)
    bhat1 /= 1 - (2 * m2 - 5) / 30 * x2

    ahat2 = 1im * x2 * (m2 - 1) / 15 * (1 - x2 / 14)
    ahat2 /= 2 * m2 + 3 - (2 * m2 - 7) / 14 * x2

    S1 = @. 1.5 * x3 * (ahat1 + bhat1 * μ + 5 / 3 * ahat2 * μ)
    S2 = @. 1.5 * x3 * (bhat1 + ahat1 * μ + 5 / 3 * ahat2 * (2 * μ^2 - 1))

    norm = √(π * 6 * x3 * real(ahat1 + bhat1 + 5 * ahat2 / 3))
    S1 ./= norm
    S2 ./= norm

    return S1, S2
end

"""
    normalization_factor(a, b, x; norm)

Figure out scattering function normalization.

# Parameters
- `m`: complex index of refraction of sphere
- `x`: dimensionless sphere size
- `norm`: symbol describing type of normalization

# Output
scaling factor needed for scattering function
"""
function normalization_factor(m, x; norm = :albedo)
    norm === :bohren && return 1 / 2
    norm === :wiscombe && return 1.0

    if norm ∈ (:qsca, :scattering_efficiency)
        return x * √(π)
    end

    qext, qsca, _, _ = mie(m, x)

    if norm ∈ (:a, :albedo)
        return x * √(π * qext)
    end

    if norm ∈ (:one, :unity)
        return x * √(π * qsca)
    end

    norm === :four_pi && return x * √(qsca / 4)

    if norm ∈ (:qext, :extinction_efficiency)
        return x * √(qsca * π / qext)
    end

    throw(
        ArgumentError(
            "Normalization must be one of :albedo (default), :one, :four_pi, :qext,
            :qsca, :bohren or :wiscombe",
        ),
    )
end

"""
Calculate the scattering amplitude functions for spheres.
The amplitude functions have been normalized so that when integrated
over all `4*π` solid angles, the integral will be `qext*pi*x^2`.
The units are weird, ``sr^{-0.5}``.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(``θ``), to calculate scattering amplitudes
- `norm` (optional): The normalization. Must be one of :albedo (default), :one, :four_pi, :qext,
:qsca, :bohren or :wiscombe"
- `use_threads` (optional): Flag whether to use threads (default: true)

# Output
`S1`, `S2`: the scattering amplitudes at each angle µ [``sr^{-0.5}``]
"""
function mie_S1_S2(m, x, μ; norm = :albedo, use_threads = true)
    a, b = mie_An_Bn(m, x)

    nangles = length(μ)
    S1 = Vector{ComplexF64}(undef, nangles)
    S2 = Vector{ComplexF64}(undef, nangles)

    nstop = length(a)
    normalization = normalization_factor(m, x; norm)
    @batch minbatch = (use_threads ? 1 : typemax(Int)) for k = 1:nangles
        pi_nm2 = 0.0
        pi_nm1 = 1.0
        @inbounds S1[k], S2[k] = zero(eltype(S1)), zero(eltype(S2))
        @inbounds for n = 1:nstop
            τ_nm1 = n * μ[k] * pi_nm1 - (n + 1) * pi_nm2
            S1[k] += (2 * n + 1) * (pi_nm1 * a[n] + τ_nm1 * b[n]) / (n + 1) / n
            S2[k] += (2 * n + 1) * (τ_nm1 * a[n] + pi_nm1 * b[n]) / (n + 1) / n
            temp = pi_nm1
            pi_nm1 = ((2 * n + 1) * μ[k] * pi_nm1 - (n + 1) * pi_nm2) / n
            pi_nm2 = temp
        end
        @inbounds S1[k] /= normalization
        @inbounds S2[k] /= normalization
    end

    return S1, S2
end

"""
    mie_cdf(m, x, num; norm = :albedo)

Create a CDF for unpolarized scattering uniformly spaced in cos(θ).
The CDF covers scattered (exit) angles ranging from 180 to 0 degrees.
(The cosines are uniformly distributed over -1 to 1.) Because the angles
are uniformly distributed in cos(theta), the scattering function is not
sampled uniformly and therefore huge array sizes are needed to adequately
sample highly anisotropic phase functions.
Since this is a cumulative distribution function, the maximum value
should be 1.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `num`: length of desired CDF array

# Output
- `µ`: array of cosines of angles
- `cdf`: array of cumulative distribution function values
"""
function mie_cdf(m, x, num; norm = :albedo)
    μ = LinRange(-1, 1, num)
    s1, s2 = mie_S1_S2(m, x, μ; norm)

    s = (abs2.(s1) + abs2.(s2)) / 2

    cdf = zeros(num)
    total = 0.0
    @inbounds for i = 1:num
        total += s[i] * 2 * π * 2 / num
        cdf[i] = total
    end

    return μ, cdf
end

"""
    mie_cdf(m, x, num; norm = :albedo)

Create a CDF for unpolarized scattering for uniform CDF.
The CDF covers scattered (exit) angles ranging from 180 to 0 degrees.
(The cosines are uniformly distributed over -1 to 1.) These angles mu
correspond to uniform spacing of the cumulative distribution function
for unpolarized Mie scattering where cdf[i] = i/(num-1).
This is a brute force implementation that solves the problem by
calculating the CDF at many points and then scanning to find the
specific angles that correspond to uniform interval of the CDF.
Since this is a cumulative distribution function, the maximum value
should be 1.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `num`: length of desired CDF array

# Output
- `µ`: array of cosines of angles (irregularly spaced)
- `cdf`: array of cumulative distribution function values
"""
function mie_mu_with_uniform_cdf(m, x, num; norm = :albedo)
    big_num = 2000 # large to wirk with x up to 10
    big_μ, big_cdf = mie_cdf(m, x, big_num; norm)
    μ, cdf = Vector{Float64}(undef, num), Vector{Float64}(undef, num)

    μ[1] = -1
    cdf[1] = 0

    big_k = 0
    @inbounds for k = 1:(num-2)
        target = k / (num - 1)
        while big_cdf[big_k+1] < target
            big_k += 1
        end

        Δ = big_cdf[big_k+1] - target
        Δ_cdf = big_cdf[big_k+1] - big_cdf[big_k]
        Δ_μ = big_μ[big_k+1] - big_μ[big_k]

        μ[k+1] = big_μ[big_k+1] - Δ / Δ_cdf * Δ_μ
        cdf[k+1] = target
    end

    μ[end] = 1
    cdf[end] = 1

    return μ, cdf
end

"""
    generate_mie_costheta(μ_cdf)

Generate a new scattering angle using a cdf.
A uniformly spaced cumulative distribution function (CDF) is needed.
New random angles are generated by selecting a random interval
μ[i] to μ[i+1] and choosing an angle uniformly distributed over
the interval.

# Parameters:
- `μ_cdf`: a cumulative distribution function

# Output:
-  an array of random scattering angle cosines based on the CDF supplied.
"""
function generate_mie_costheta(μ_cdf)
    index = rand(1:(length(μ_cdf)-1))

    x = μ_cdf[index]
    x += (μ_cdf[index+1] - μ_cdf[index]) * rand()

    return x
end

"""
    i_per(m, x, μ; norm = :albedo)

Return the scattered intensity in a plane normal to the incident light.
This is the scattered intensity in a plane that is perpendicular to the
field of the incident plane wave. The intensity is normalized such
that the integral of the unpolarized intensity over 4π steradians
is equal to the single scattering albedo.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(theta), to calculate intensities

# Output
The intensity at each angle in the array µ.  Units [1/sr]
"""
function i_per(m, x, μ; norm = :albedo)
    s1, _ = mie_S1_S2(m, x, μ; norm)
    return abs2.(s1)
end

"""
    i_par(m, x, μ; norm = :albedo)

Return the scattered intensity in a plane parallel to the incident light.
This is the scattered intensity in a plane that is perpendicular to the
field of the incident plane wave. The intensity is normalized such
that the integral of the unpolarized intensity over 4π steradians
is equal to the single scattering albedo.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(theta), to calculate intensities

# Output
The intensity at each angle in the array µ.  Units [1/sr]
"""
function i_par(m, x, μ; norm = :albedo)
    _, s2 = mie_S1_S2(m, x, μ; norm)
    return abs2.(s2)
end

"""
    i_unpolarized(m, x, μ; norm = :albedo)

Return the unpolarized scattered intensity at specified angles.
This is the average value for randomly polarized incident light.
The intensity is normalized such
that the integral of the unpolarized intensity over 4π steradians
is equal to the single scattering albedo.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `µ`: the angles, cos(theta), to calculate intensities

# Output
The intensity at each angle in the array µ.  Units [1/sr]
    """
function i_unpolarized(m, x, μ; norm = :albedo)
    s1, s2 = mie_S1_S2(m, x, μ; norm)
    return @. (abs2(s1) + abs2(s2)) / 2
end

"""
    ez_mie(m, d, λ0, n_env = 1.0)

Calculate the efficiencies of a sphere.

# Parameters
- `m`: the complex index of refraction of the sphere    [-]
- `d`: the diameter of the sphere                       [same units as lambda0]
- `λ0`: wavelength in a vacuum                          [same units as d]
- `n_env`: real index of medium around sphere, optional.
- `use_threads` (optional): Flag whether to use threads (default: true)

# Output
- `qext`: the total extinction efficiency                  [-]
- `qsca`: the scattering efficiency                        [-]
- `qback`: the backscatter efficiency                      [-]
- `g`: the average cosine of the scattering phase function [-]
"""
function ez_mie(m, d, λ0, n_env = 1.0; use_threads = true)
    m_env = @. m / n_env
    x_env = @. π * d / (λ0 / n_env)
    return mie(m_env, x_env; use_threads)
end

"""
    ez_intensities(m, d, λ0, μ, n_env = 1.0, norm = :albedo)

Return the scattered intensities from a sphere.
These are the scattered intensities in a plane that is parallel (ipar) and
perpendicular (iper) to the field of the incident plane wave.
The scattered intensity is normalized such that the integral of the
unpolarized intensity over 4𝜋 steradians is equal to the single scattering
albedo.  The scattered intensity has units of inverse steradians [1/sr].
The unpolarized scattering is the average of the two scattered intensities.

# Parameters
- `m`: the complex index of refraction of the sphere    [-]
- `d`: the diameter of the sphere                       [same units as lambda0]
- `λ0`: wavelength in a vacuum                          [same units as d]
- `µ`: the cos(θ) of each direction desired             [-]
- `n_env`: real index of medium around sphere, optional.

# Output
`ipar`, `iper`: scattered intensity in parallel and perpendicular planes [1/sr]
"""
function ez_intensities(m, d, λ0, μ, n_env = 1.0, norm = :albedo)
    m_env = m / n_env
    λ_env = λ0 / n_env
    x_env = π * d / λ_env
    s1, s2 = mie_S1_S2(m_env, x_env, μ; norm)
    ipar = abs2(s2)
    iper = abs2(s1)
    return ipar, iper
end

"""
    mie_phase_matrix(m, x, μ; norm=:albedo)

Calculate the phase scattering matrix.

The units are ``sr^{-1}``.
The phase scattering matrix is computed from the scattering amplitude
functions, according to equations 5.2.105-6 in K. N. Liou (**2002**) -
*An Introduction to Atmospheric Radiation*, Second Edition.

# Parameters
- `m`: the complex index of refraction of the sphere
- `x`: the size parameter of the sphere
- `μ`: the angles, cos(theta), at which to calculate the phase scattering matrix

# Output
- `p`: The phase scattering matrix [``sr^{-1}``]
"""
function mie_phase_matrix end

function mie_phase_matrix(m, x, μ::AbstractVector; norm = :albedo)
    s1, s2 = mie_S1_S2(m, x, μ; norm)

    phase = zeros(4, 4, length(μ))
    for i in eachindex(s1)
        s1_star = conj(s1[i])
        s2_star = conj(s2[i])

        m1 = abs2(s1[i])
        m2 = abs2(s2[i])
        s21 = real(0.5 * (s1[i] * s2_star + s2[i] * s1_star))
        d21 = real(-0.5*im * (s1[i] * s2_star - s2[i] * s1_star))

        phase[1, 1, i] = 0.5 * (m2 + m1)
        phase[1, 2, i] = 0.5 * (m2 - m1)
        phase[2, 1, i] = phase[1, 2, i]
        phase[2, 2, i] = phase[1, 1, i]
        phase[3, 3, i] = s21
        phase[3, 4, i] = -d21
        phase[4, 3, i] = d21
        phase[4, 4, i] = s21
    end

    return phase
end

mie_phase_matrix(m, x, μ::Number; norm = :albedo) =
    reshape(mie_phase_matrix(m, x, [µ]; norm), 4, 4)
