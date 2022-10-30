
function mie end

function mie(m::T, x::V) where {T <: Number, V <: Number}
    TT = float(real(promote_type(T, V)))
    if real(m) == zero(real(T)) && x < TT(0.1)
        qext, qsca, qback, g = small_conducting_mie(m, x)
    elseif real(m) > zero(real(T)) && abs(m) * x < TT(0.1)
        qext, qsca, qback, g = small_mie(m, x)
    else
        a, b = mie_An_Bn(m, x)

        nmax = length(a)
        n = 1:nmax
        cn = 2.0 * n .+ 1.0

        qext = 2 * sum(@.(cn * (real(a) + real(b)))) / x^2
        qsca = qext

        if imag(m) != zero(TT)
            qsca = 2 * sum(@.(cn * (abs2(a) + abs2(b)))) / x^2
        end

        qback = abs(sum(@.((-1)^n * cn * (a - b))))^2 / x^2

        c1n = @. n * (n + 2) / (n + 1)
        c2n = @. cn / n / (n + 1)
        g = zero(TT)
        for i=1:(nmax-1)
            asy1 = c1n[i] * real(a[i] * conj(a[i+1]) + b[i] * conj(b[i+1]))
            asy2 = c2n[i] * real(a[i] * conj(b[i]))
            g += 4 * (asy1 + asy2) / qsca / x^2
        end
    end

    return qext, qsca, qback, g
end

function mie(m::AbstractVector{T}, x::AbstractVector{V}, force_threads=false) where {T <: Number, V <: Number}
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

    if force_threads || len > 50
        Threads.@threads for i in eachindex(g)
            mm = mlen > 1 ? m[i] : first(m)
            xx = xlen > 1 ? x[i] : first(x)
            qext[i], qsca[i], qback[i], g[i] = mie(mm, xx)
        end
    else
        for i in eachindex(g)
            mm = mlen > 1 ? m[i] : first(m)
            xx = xlen > 1 ? x[i] : first(x)
            qext[i], qsca[i], qback[i], g[i] = mie(mm, xx)
        end
    end

    return qext, qsca, qback, g
end
mie(m::Number, x::AbstractVector) = mie([m], x)
mie(m::AbstractVector, x::Number) = mie(m, [x])

function small_mie(m, x)
    m2 = m^2
    x2 = x^2

    D = complex(m2 + 2 + (1 - 0.7 * m2) * x2)
    D -= (8 * m^4 - 385*m2 + 350) * x^4 / 1400.0
    D += 2*im * (m2 - 1) * x^3 * (1-0.1*x2)/3
    ahat1 = 2im * (m2 - 1) / 3 * (1-0.1*x2+(4*m2+5)*x^4/1400)/D
    bhat1 = 1im * x2 * (m2 - 1) / 45 * (1 + (2*m2 - 5) / 70 * x2)
    bhat1 /= 1 - (2*m2 - 5) / 30 * x2
    ahat2 = im *x2 * (m2 - 1)/ 15 * (1-x2/14)
    ahat2 /= 2*m2 + 3 - (2*m2 - 7)/14*x2

    T = abs2(ahat1) + abs2(bhat1) + 5/3 * abs2(ahat2)
    temp = ahat2 + bhat1
    g = real(ahat1 * conj(temp)) / T

    qsca = 6*x^4*T

    if real(m) == 0
        qext = qsca
    else
        qext = 6*x*real(ahat1 + bhat1 + 5*ahat2/3)
    end

    sback = 1.5*x^3*(ahat1 - bhat1 - 5*ahat2/3)
    qback = 4*abs2(sback)/x2

    return qext, qsca, qback, g
end

function small_conducting_mie(m, x)
    ahat1 = (im*2/3*(1-0.2*x^2)) / (1-0.5*x^2 + im*2/3*x^3)
    bhat1 = (im*(x^2-10.0)/30.0) / (1+0.5*x^2 - im*x^3/3.0)
    ahat2 = im*x^2/30
    bhat2 = -im*x^2/45

    qsca = x^4 * (6*abs2(ahat1) + 6*abs2(bhat1) + 10*abs2(ahat2) + 10*abs2(bhat2))
    qext = qsca
    g = imag(ahat1) * (imag(ahat2) + imag(bhat1))
    g += imag(bhat2) * (5/9 * imag(ahat2) + imag(bhat1))
    g += real(ahat1) * real(bhat1)
    g *= 6*x^4/qsca

    qback = 9*x^4*abs(ahat1-bhat1-5/3*(ahat2-bhat2))^2

    return qext, qsca, qback, g
end

function mie_An_Bn(m, x)
    nstop = floor(Int, x + 4.05*x^0.33333 + 2.0) + 1
    a = zeros(ComplexF64, nstop - 1)
    b = zeros(ComplexF64, nstop - 1)

    psi_nm1 = sin(x)
    psi_n = psi_nm1 / x - cos(x)
    xi_nm1 = psi_nm1 + im*cos(x)
    xi_n = psi_n + im * (cos(x)/x + sin(x))

    if real(m) > 0.0
        D = D_calc(m, x, nstop + 1)

        for n=1:(nstop-1)
            temp = D[n] / m + n / x
            a[n] = (temp * psi_n - psi_nm1) / (temp * xi_n - xi_nm1)
            temp = D[n] * m + n / x
            b[n] = (temp * psi_n - psi_nm1) / (temp * xi_n - xi_nm1)
            xi = (2*n + 1) * xi_n / x - xi_nm1
            xi_nm1 = xi_n
            xi_n = xi
            psi_nm1 = psi_n
            psi_n = real(xi_n)
        end
    else
        for n=1:(nstop-1)
            a[n] = (n*psi_n/x - psi_nm1) /(n*xi_n/x-xi_nm1)
            b[n] = psi_n/xi_n
            xi = (2*n+1)*xi_n/x - xi_nm1
            xi_nm1 = xi_n
            xi_n = xi
            psi_nm1 = psi_n
            psi_n = real(xi_n)
        end
    end

    return a, b
end

function D_calc(m, x, N)
    n = real(m)
    κ = abs(imag(m))
    D = zeros(ComplexF64, N)

    if n < 1 || n > 10 || κ > 10 || (x * κ) >= (3.9 - 10.8*n + 13.78*n^2)
        D_downwards!(m*x, N, D)
    else
        D_upwards!(m*x, N, D)
    end

    return D
end

function D_upwards!(z, N, D)
    _exp = exp(-2im*z)
    D[1] = -1/z+(1-_exp) / ((1-_exp)/z-im*(1+_exp))
    for n=2:N
        D[n] = 1/(n/z-D[n-1])-n/z
    end
end

function D_downwards!(z, N, D)
    last_D = Lentz_Dn(z, N)
    for n=N:-1:2
        last_D = n/z-1.0/(last_D+n/z)
        D[n-1] = last_D
    end
end

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

function small_mie_conducting_S1_S2(m, x, μ)
    ahat1 = 2im/3*(1-0.2*x^2)/(1-0.5*x^2+2im/3*x^3)
    bhat1 = 1im/3*(0.1*x^2-1)/(1+0.5*x^2-1im/3*x^3)
    ahat2 = 1im/30*x^2
    bhat2 = -1im*x^2/45

    S1 = @. 1.5*x^3*(ahat1+bhat1*μ+5/3*ahat2*μ+5/3*bhat2*(2*μ^2-1))
    S2 = @. 1.5*x^3*(bhat1+ahat1*μ+5/3*bhat2*μ+5/3*ahat2*(2*μ^2-1))

    qext = x^4*(6*abs2(ahat1) + 6 * abs(bhat1) + 10*abs2(ahat2) + 10*abs2/bhat2)

    norm = √(qext*π*x^2)
    S1 /= norm
    S2 /= norm

    return S1, S2
end
