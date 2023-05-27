var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [MieScattering]","category":"page"},{"location":"api/#MieScattering.MieScattering","page":"API","title":"MieScattering.MieScattering","text":"Mie scattering calculations for perfect spheres based on miepython. Extensive documentation is at (https://miepython.readthedocs.io).\n\nMieScattering.jl is a Julia package to calculate light scattering of a plane wave by non-np.absorbing, partially-np.absorbing, or perfectly conducting spheres.\n\nThe extinction efficiency, scattering efficiency, backscattering, and scattering asymmetry for a sphere with complex index of refraction m, diameter d, and wavelength lambda can be found by:\n\n    qext, qsca, qback, g = ez_mie(m, d, λ0)\n\nThe normalized scattering values for angles µ=cos(θ) are:\n\n    Ipar, Iper = ez_intensities(m, d, λ0, µ)\n\nIf the size parameter is known, then use:\n\n    mie(m, x)\n\nMie scattering amplitudes S1 and S2 (complex numbers):\n\n    mie_S1_S2(m, x, μ)\n\nNormalized Mie scattering intensities for angles µ=cos(θ):\n\n    i_per(m, x, µ)\n    i_par(m, x, µ)\n    i_unpolarized(m, x, µ)\n\n\n\n\n\n","category":"module"},{"location":"api/#MieScattering.D_calc-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.D_calc","text":"D_calc(m, x, N)\n\nCompute the logarithmic derivative using best method.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nN: order of Ricatti-Bessel function\n\nOutput\n\nThe values of the Ricatti-Bessel function for orders from 0 to N.\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.D_downwards!-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.D_downwards!","text":"D_downwards!(z, N, D)\n\nCompute the logarithmic derivative by downwards recurrence.\n\nParameters\n\nz: function argument\nN: order of Ricatti-Bessel function\nD: gets filled with the Ricatti-Bessel function values for orders      from 0 to N for an argument z using the downwards recurrence relations.\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.D_upwards!-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.D_upwards!","text":"D_upwards!(z, N, D)\n\nCompute the logarithmic derivative by upwards recurrence.\n\nParameters\n\nz: function argument\nN: order of Ricatti-Bessel function\nD: gets filled with the Ricatti-Bessel function values for orders      from 0 to N for an argument z using the upwards recurrence relations.\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.Lentz_Dn-Tuple{Any, Any}","page":"API","title":"MieScattering.Lentz_Dn","text":"Lentz_Dn(z, N)\n\nCompute the logarithmic derivative of the Ricatti-Bessel function.\n\nParameters\n\nz: function argument\nN: order of Ricatti-Bessel function\n\nOutput\n\nThis returns the Ricatti-Bessel function of order N with argument z using the continued fraction technique of Lentz, Appl. Opt., 15, 668-671, (1976).\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.ez_intensities","page":"API","title":"MieScattering.ez_intensities","text":"ez_intensities(m, d, λ0, μ, n_env = 1.0, norm = :albedo)\n\nReturn the scattered intensities from a sphere. These are the scattered intensities in a plane that is parallel (ipar) and perpendicular (iper) to the field of the incident plane wave. The scattered intensity is normalized such that the integral of the unpolarized intensity over 4𝜋 steradians is equal to the single scattering albedo.  The scattered intensity has units of inverse steradians [1/sr]. The unpolarized scattering is the average of the two scattered intensities.\n\nParameters\n\nm: the complex index of refraction of the sphere    [-]\nd: the diameter of the sphere                       [same units as lambda0]\nλ0: wavelength in a vacuum                          [same units as d]\nµ: the cos(θ) of each direction desired             [-]\nn_env: real index of medium around sphere, optional.\n\nOutput\n\nipar, iper: scattered intensity in parallel and perpendicular planes [1/sr]\n\n\n\n\n\n","category":"function"},{"location":"api/#MieScattering.ez_mie","page":"API","title":"MieScattering.ez_mie","text":"ez_mie(m, d, λ0, n_env = 1.0)\n\nCalculate the efficiencies of a sphere.\n\nParameters\n\nm: the complex index of refraction of the sphere    [-]\nd: the diameter of the sphere                       [same units as lambda0]\nλ0: wavelength in a vacuum                          [same units as d]\nn_env: real index of medium around sphere, optional.\nuse_threads (optional): Flag whether to use threads (default: true)\n\nOutput\n\nqext: the total extinction efficiency                  [-]\nqsca: the scattering efficiency                        [-]\nqback: the backscatter efficiency                      [-]\ng: the average cosine of the scattering phase function [-]\n\n\n\n\n\n","category":"function"},{"location":"api/#MieScattering.generate_mie_costheta-Tuple{Any}","page":"API","title":"MieScattering.generate_mie_costheta","text":"generate_mie_costheta(μ_cdf)\n\nGenerate a new scattering angle using a cdf. A uniformly spaced cumulative distribution function (CDF) is needed. New random angles are generated by selecting a random interval μ[i] to μ[i+1] and choosing an angle uniformly distributed over the interval.\n\nParameters:\n\nμ_cdf: a cumulative distribution function\n\nOutput:\n\nan array of random scattering angle cosines based on the CDF supplied.\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.i_par-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.i_par","text":"i_par(m, x, μ; norm = :albedo)\n\nReturn the scattered intensity in a plane parallel to the incident light. This is the scattered intensity in a plane that is perpendicular to the field of the incident plane wave. The intensity is normalized such that the integral of the unpolarized intensity over 4π steradians is equal to the single scattering albedo.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(theta), to calculate intensities\n\nOutput\n\nThe intensity at each angle in the array µ.  Units [1/sr]\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.i_per-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.i_per","text":"i_per(m, x, μ; norm = :albedo)\n\nReturn the scattered intensity in a plane normal to the incident light. This is the scattered intensity in a plane that is perpendicular to the field of the incident plane wave. The intensity is normalized such that the integral of the unpolarized intensity over 4π steradians is equal to the single scattering albedo.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(theta), to calculate intensities\n\nOutput\n\nThe intensity at each angle in the array µ.  Units [1/sr]\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.i_unpolarized-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.i_unpolarized","text":"i_unpolarized(m, x, μ; norm = :albedo)\n\nReturn the unpolarized scattered intensity at specified angles. This is the average value for randomly polarized incident light. The intensity is normalized such that the integral of the unpolarized intensity over 4π steradians is equal to the single scattering albedo.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(theta), to calculate intensities\n\nOutput\n\nThe intensity at each angle in the array µ.  Units [1/sr]\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.mie","page":"API","title":"MieScattering.mie","text":"mie(m, x)\n\nCalculate the efficiencies for a sphere where m or x may be vectors.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\n\nOutput\n\nqext: the total extinction efficiency\nqsca: the scattering efficiency\nqback: the backscatter efficiency\ng: the average cosine of the scattering phase function\n\n\n\n\n\n","category":"function"},{"location":"api/#MieScattering.mie_An_Bn-Tuple{Any, Any}","page":"API","title":"MieScattering.mie_An_Bn","text":"mie_An_Bn(m, x)\n\nCompute arrays of Mie coefficients A and B for a sphere. This estimates the size of the arrays based on Wiscombe's formula. The length of the arrays is chosen so that the error when the series are summed is around 1e-6.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\n\nOutput\n\nAn, Bn: arrays of Mie coefficents\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.mie_S1_S2-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.mie_S1_S2","text":"Calculate the scattering amplitude functions for spheres. The amplitude functions have been normalized so that when integrated over all 4*π solid angles, the integral will be qext*pi*x^2. The units are weird, sr^-05.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(θ), to calculate scattering amplitudes\nnorm (optional): The normalization. Must be one of :albedo (default), :one, :four_pi, :qext,\n\n:qsca, :bohren or :wiscombe\"\n\nuse_threads (optional): Flag whether to use threads (default: true)\n\nOutput\n\nS1, S2: the scattering amplitudes at each angle µ [sr^-05]\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.mie_cdf-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.mie_cdf","text":"mie_cdf(m, x, num; norm = :albedo)\n\nCreate a CDF for unpolarized scattering uniformly spaced in cos(θ). The CDF covers scattered (exit) angles ranging from 180 to 0 degrees. (The cosines are uniformly distributed over -1 to 1.) Because the angles are uniformly distributed in cos(theta), the scattering function is not sampled uniformly and therefore huge array sizes are needed to adequately sample highly anisotropic phase functions. Since this is a cumulative distribution function, the maximum value should be 1.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nnum: length of desired CDF array\n\nOutput\n\nµ: array of cosines of angles\ncdf: array of cumulative distribution function values\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.mie_mu_with_uniform_cdf-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.mie_mu_with_uniform_cdf","text":"mie_cdf(m, x, num; norm = :albedo)\n\nCreate a CDF for unpolarized scattering for uniform CDF. The CDF covers scattered (exit) angles ranging from 180 to 0 degrees. (The cosines are uniformly distributed over -1 to 1.) These angles mu correspond to uniform spacing of the cumulative distribution function for unpolarized Mie scattering where cdf[i] = i/(num-1). This is a brute force implementation that solves the problem by calculating the CDF at many points and then scanning to find the specific angles that correspond to uniform interval of the CDF. Since this is a cumulative distribution function, the maximum value should be 1.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nnum: length of desired CDF array\n\nOutput\n\nµ: array of cosines of angles (irregularly spaced)\ncdf: array of cumulative distribution function values\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.mie_phase_matrix","page":"API","title":"MieScattering.mie_phase_matrix","text":"mie_phase_matrix(m, x, μ; norm=:albedo)\n\nCalculate the phase scattering matrix.\n\nThe units are sr^-1. The phase scattering matrix is computed from the scattering amplitude functions, according to equations 5.2.105-6 in K. N. Liou (2002) - An Introduction to Atmospheric Radiation, Second Edition.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nμ: the angles, cos(theta), at which to calculate the phase scattering matrix\n\nOutput\n\np: The phase scattering matrix [sr^-1]\n\n\n\n\n\n","category":"function"},{"location":"api/#MieScattering.normalization_factor-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.normalization_factor","text":"normalization_factor(a, b, x; norm)\n\nFigure out scattering function normalization.\n\nParameters\n\na: complex array of An coefficients\nb: complex array of Bn coefficients\nx: dimensionless sphere size\nnorm: symbol describing type of normalization\n\nOutput\n\nscaling factor needed for scattering function\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.small_conducting_mie-Tuple{Any}","page":"API","title":"MieScattering.small_conducting_mie","text":"small_conducting_mie(m, x)\n\nCalculate the efficiencies for a small conducting spheres. Typically used for small conducting spheres where x < 0.1 and real(m) == 0.\n\nParameters\n\nx: the size parameter of the sphere\n\nOutput\n\nqext: the total extinction efficiency\nqsca: the scattering efficiency\nqback: the backscatter efficiency\ng: the average cosine of the scattering phase function\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.small_mie-Tuple{Any, Any}","page":"API","title":"MieScattering.small_mie","text":"small_mie(m, x)\n\nCalculate the efficiencies for a small sphere. Typically used for small spheres where x<0.1\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\n\nOutput\n\nqext: the total extinction efficiency\nqsca: the scattering efficiency\nqback: the backscatter efficiency\ng: the average cosine of the scattering phase function\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.small_mie_S1_S2-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.small_mie_S1_S2","text":"Calculate the scattering amplitude functions for small spheres (x<0.1). The amplitude functions have been normalized so that when integrated over all 4*π solid angles, the integral will be qext*pi*x^2. The units are weird, sr^-05\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(θ), to calculate scattering amplitudes\n\nOutput\n\nS1, S2: the scattering amplitudes at each angle µ [sr^-05]\n\n\n\n\n\n","category":"method"},{"location":"api/#MieScattering.small_mie_conducting_S1_S2-Tuple{Any, Any, Any}","page":"API","title":"MieScattering.small_mie_conducting_S1_S2","text":"small_mie_conducting_S1_S2(m, x, μ)\n\nCalculate the scattering amplitudes for small conducting spheres. The spheres are small perfectly conducting (reflecting) spheres (x<0.1). The amplitude functions have been normalized so that when integrated over all 4𝜋 solid angles, the integral will be qext(𝜋x²). The units are weird, sr^-05.\n\nParameters\n\nm: the complex index of refraction of the sphere\nx: the size parameter of the sphere\nµ: the angles, cos(θ), to calculate scattering amplitudes\n\nOutput\n\nS1, S2: the scattering amplitudes at each angle µ [sr^-05]\n\n\n\n\n\n","category":"method"},{"location":"performance/#Performance-in-comparison-to-miepython","page":"Performance","title":"Performance in comparison to miepython","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"The following snippets show how MieScattering.jl performs compared to miepython. The tests are performed by calling miepython using PyCall. This might seem as an unfair comparison due to the (little) overhead of calling Python from Julia. However, if you are reading this you already want to use Julia for this task so resorting to PyCall might be your way to go anyway if there is no native Julia solution.","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"The performance is compared vs the jit-ed and non-jit-ed version. Additionally the Julia version is used with the parameter use_threads=false  in order to show the difference when using multi-threading in Julia (using the great package Polyester.jl). Multi-threading is the default setting, even though there is a slight disadvantage for small workloads. Even in single-threaded operation MieScattering.jl gives slightly better performance compared to the numba-accelerated miepython version (which essentially also does LLVM compilation of the Python code).","category":"page"},{"location":"performance/#Version-info","page":"Performance","title":"Version info","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"Julia Version 1.8.2\nCommit 36034abf260 (2022-09-29 15:21 UTC)\nPlatform Info:\n  OS: macOS (arm64-apple-darwin21.3.0)\n  CPU: 10 × Apple M1 Pro\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-13.0.1 (ORCJIT, apple-m1)\n  Threads: 8 on 8 virtual cores\nEnvironment:\n  JULIA_EDITOR = code\n\nsys.version = \"3.10.6 | packaged by conda-forge | (main, Aug 22 2022, 20:40:44) [Clang 13.0.1 ]\"","category":"page"},{"location":"performance/#Size-parameter","page":"Performance","title":"Size parameter","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"using MieScattering\nusing BenchmarkTools\nusing PyCall\nusing Conda\nConda.pip_interop(true)\nConda.pip(\"install\", \"miepython\")\nusing CairoMakie\n\nmiepython_jit = pyimport(\"miepython.miepython\")\nmiepython = pyimport(\"miepython.miepython_nojit\")\n\nntests = 6\nm = 1.5\nN = round.(Int, 10 .^ (LinRange(0, 3, ntests)))\nN[1] = 2\n\nresult = zeros(ntests)\nresult_threads = zeros(ntests)\nresult_python = zeros(ntests)\nresult_python_jit = zeros(ntests)\n\nfor i = 1:ntests\n    x = collect(LinRange(0.1, 20.0, N[i]))\n    a = @benchmark mie($m, $x, use_threads=false)\n    result[i] = median(a).time\n    a = @benchmark mie($m, $x)\n    result_threads[i] = median(a).time\n    a = @benchmark miepython.mie($m, $x)\n    result_python[i] = median(a).time\n    a = @benchmark miepython_jit.mie($m, $x)\n    result_python_jit[i] = median(a).time\nend\n\nimp = result_python_jit ./ result_threads\n\nminimp, maximp = extrema(imp)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        yscale = log10,\n        xscale = log10,\n        xlabel = \"Number of sphere sizes calculated\",\n        ylabel = \"Execution Time [s]\",\n        title = \"Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)\",\n    )\nscatterlines!(N, result * 1e-9, label = \"Julia single-threaded (MieScattering.jl)\")\nscatterlines!(N, result_threads * 1e-9, label = \"Julia multi-threaded (MieScattering.jl)\")\nscatterlines!(N, result_python * 1e-9, label = \"PyCall (miepython)\")\nscatterlines!(N, result_python_jit * 1e-9, label = \"PyCall (miepython-jit)\")\naxislegend(ax; position = :lt)\nfig","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"(Image: )","category":"page"},{"location":"performance/#Embedded-spheres","page":"Performance","title":"Embedded spheres","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"ntests = 6\nmwater = 4/3\nm = 1.0\nmm = m/mwater\nr = 500 # nm\n\nN = round.(Int, 10 .^ (LinRange(0, 3, ntests)))\nN[1] = 2\n\nresult = zeros(ntests)\nresult_threads = zeros(ntests)\nresult_python = zeros(ntests)\nresult_python_jit = zeros(ntests)\n\nfor i = 1:ntests\n    λ0 = collect(LinRange(300, 800, N[i]))\n    xx = 2π*r*mwater./λ0\n    a = @benchmark mie($mm, $xx, use_threads=false)\n    result[i] = median(a).time\n    a = @benchmark mie($mm, $xx)\n    result_threads[i] = median(a).time\n    a = @benchmark miepython.mie($mm, $xx)\n    result_python[i] = median(a).time\n    a = @benchmark miepython_jit.mie($mm, $xx)\n    result_python_jit[i] = median(a).time\nend\n\nimp = result_python_jit ./ result_threads\n\nminimp, maximp = extrema(imp)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        yscale = log10,\n        xscale = log10,\n        xlabel = \"Number of wavelengths calculated\",\n        ylabel = \"Execution Time [s]\",\n        title = \"Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)\",\n    )\nscatterlines!(N, result * 1e-9, label = \"Julia single-threaded (MieScattering.jl)\")\nscatterlines!(N, result_threads * 1e-9, label = \"Julia multi-threaded (MieScattering.jl)\")\nscatterlines!(N, result_python * 1e-9, label = \"PyCall (miepython)\")\nscatterlines!(N, result_python_jit * 1e-9, label = \"PyCall (miepython-jit)\")\naxislegend(ax; position = :lt)\nfig","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"(Image: )","category":"page"},{"location":"performance/#ez_mie-testing","page":"Performance","title":"ez_mie testing","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"ntests = 6\nm_sphere = 1.0\nn_water = 4/3\nd = 1000 # nm\n\nN = round.(Int, 10 .^ (LinRange(0, 3, ntests)))\nN[1] = 2\n\nresult = zeros(ntests)\nresult_threads = zeros(ntests)\nresult_python = zeros(ntests)\nresult_python_jit = zeros(ntests)\n\nfor i = 1:ntests\n    λ0 = collect(LinRange(300, 800, N[i]))\n    a = @benchmark ez_mie($m_sphere, $d, $λ0, $n_water, use_threads=false)\n    result[i] = median(a).time\n    a = @benchmark ez_mie($m_sphere, $d, $λ0, $n_water)\n    result_threads[i] = median(a).time\n    a = @benchmark miepython.ez_mie($m_sphere, $d, $λ0, $n_water)\n    result_python[i] = median(a).time\n    a = @benchmark miepython_jit.ez_mie($m_sphere, $d, $λ0, $n_water)\n    result_python_jit[i] = median(a).time\nend\n\nimp = result_python_jit ./ result_threads\n\nminimp, maximp = extrema(imp)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        yscale = log10,\n        xscale = log10,\n        xlabel = \"Number of wavelengths calculated\",\n        ylabel = \"Execution Time [s]\",\n        title = \"Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)\",\n    )\nscatterlines!(N, result * 1e-9, label = \"Julia single-threaded (MieScattering.jl)\")\nscatterlines!(N, result_threads * 1e-9, label = \"Julia multi-threaded (MieScattering.jl)\")\nscatterlines!(N, result_python * 1e-9, label = \"PyCall (miepython)\")\nscatterlines!(N, result_python_jit * 1e-9, label = \"PyCall (miepython-jit)\")\naxislegend(ax; position = :lt)\nfig","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"(Image: )","category":"page"},{"location":"performance/#Scattering-Phase-Function","page":"Performance","title":"Scattering Phase Function","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"ntests = 6\nm = 1.5\nx = π/3\n\nN = round.(Int, 10 .^ (LinRange(0, 3, ntests)))\nN[1] = 2\n\nresult = zeros(ntests)\nresult_threads = zeros(ntests)\nresult_python = zeros(ntests)\nresult_python_jit = zeros(ntests)\n\nfor i = 1:ntests\n    θ = LinRange(-180, 180, N[i])\n    μ = cosd.(θ)\n    a = @benchmark mie_S1_S2($m, $x, $μ, use_threads=false)\n    result[i] = median(a).time\n    a = @benchmark mie_S1_S2($m, $x, $μ)\n    result_threads[i] = median(a).time\n    a = @benchmark miepython.mie_S1_S2($m, $x, $μ)\n    result_python[i] = median(a).time\n    a = @benchmark miepython_jit.mie_S1_S2($m, $x, $μ)\n    result_python_jit[i] = median(a).time\nend\n\nimp = result_python_jit ./ result_threads\n\nminimp, maximp = extrema(imp)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        yscale = log10,\n        xscale = log10,\n        xlabel = \"Number of angles calculated\",\n        ylabel = \"Execution Time [s]\",\n        title = \"Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)\",\n    )\nscatterlines!(N, result * 1e-9, label = \"Julia single-threaded (MieScattering.jl)\")\nscatterlines!(N, result_threads * 1e-9, label = \"Julia multi-threaded (MieScattering.jl)\")\nscatterlines!(N, result_python * 1e-9, label = \"PyCall (miepython)\")\nscatterlines!(N, result_python_jit * 1e-9, label = \"PyCall (miepython-jit)\")\naxislegend(ax; position = :lt)\nfig","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"(Image: )","category":"page"},{"location":"performance/#As-a-function-of-sphere-size","page":"Performance","title":"As a function of sphere size","text":"","category":"section"},{"location":"performance/","page":"Performance","title":"Performance","text":"ntests = 6\nm = 1.5 - 0.1*im\n\nx = round.(Int, 10 .^ (LinRange(0, 3, ntests)))\n\nθ = LinRange(-180, 180, 50)\nμ = cosd.(θ)\n\nresult = zeros(ntests)\nresult_threads = zeros(ntests)\nresult_python = zeros(ntests)\nresult_python_jit = zeros(ntests)\n\nfor i = 1:ntests\n    a = @benchmark mie_S1_S2($m, $x[$i], $μ, use_threads=false)\n    result[i] = median(a).time\n    a = @benchmark mie_S1_S2($m, $x[$i], $μ)\n    result_threads[i] = median(a).time\n    a = @benchmark miepython.mie_S1_S2($m, $x[$i], $μ)\n    result_python[i] = median(a).time\n    a = @benchmark miepython_jit.mie_S1_S2($m, $x[$i], $μ)\n    result_python_jit[i] = median(a).time\nend\n\nimp = result_python_jit ./ result_threads\n\nminimp, maximp = extrema(imp)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        yscale = log10,\n        xscale = log10,\n        xlabel = \"Sphere size parameter\",\n        ylabel = \"Execution Time [s]\",\n        title = \"Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)\",\n    )\nscatterlines!(N, result * 1e-9, label = \"Julia single-threaded (MieScattering.jl)\")\nscatterlines!(N, result_threads * 1e-9, label = \"Julia multi-threaded (MieScattering.jl)\")\nscatterlines!(N, result_python * 1e-9, label = \"PyCall (miepython)\")\nscatterlines!(N, result_python_jit * 1e-9, label = \"PyCall (miepython-jit)\")\naxislegend(ax; position = :lt)\nfig","category":"page"},{"location":"performance/","page":"Performance","title":"Performance","text":"(Image: )","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = MieScattering","category":"page"},{"location":"#MieScattering.jl","page":"Introduction","title":"MieScattering.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The package MieScattering.jl is a pure-Julia counterpart to the miepython package developed by Scott Prahl. This package is a 1:1 port and therefore copies the syntax, parts of the documentation and the Mie scattering related assumptions. Checkout the documentation of miepython to find out more.","category":"page"},{"location":"#Switch-from-miepython-to-MieScattering.jl","page":"Introduction","title":"Switch from miepython to MieScattering.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The syntax is mostly the same. Just convert your calls from Python to Julia and you should get the same numerical results.","category":"page"},{"location":"#Script-examples","page":"Introduction","title":"Script examples","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The following shows how to convert the four miepython example scripts to Julia.","category":"page"},{"location":"#Simple-Dielectric","page":"Introduction","title":"Simple Dielectric","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"\"\"\"\nPlot the extinction efficiency as a function of particle size.\n\nThis is a comparision of total extinction for non-absorbing and absorbing\nspheres.\n\"\"\"\n\nusing MieScattering\nusing CairoMakie\n\nx = LinRange(0.1, 100, 300)\n\n# mie() will automatically try to do the right thing\nqext, qsca, qback, g = mie(1.5, x)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        xlabel = \"Size parameter (-)\",\n        ylabel = \"Qext\",\n        title = \"Comparison of extinction for absorbing and non-absorbing spheres\",\n    )\nlines!(x, qext, color = :red, label = \"1.5\")\n\nqext, qsca, qback, g = mie(1.5 - 0.1 * im, x)\nlines!(x, qext, color = :blue, label = \"1.5-0.1j\")\n\naxislegend()\n\nfig","category":"page"},{"location":"#Glass-Spheres","page":"Introduction","title":"Glass Spheres","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"\"\"\"\nPlot the scattering efficiency for 4 micron glass spheres.\n\nThis graph shows scattering as a function of wavelength.\n\"\"\"\n\nusing MieScattering\nusing CairoMakie\n\nradius = 2 # in microns\nλ0 = LinRange(0.2, 1.2, 200)\nx = 2π * radius ./ λ0\n\n# from https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT\nm2 = @. 1 + 1.03961212 / (1 - 0.00600069867 / λ0^2)\nm2 += @. 0.231792344 / (1 - 0.0200179144 / λ0^2)\nm2 += @. 1.01046945 / (1 - 103.560653 / λ0^2)\nm = .√(m2)\nqext, qsca, qback, g = mie(m, x)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        xlabel = \"Wavelength [nm]\",\n        ylabel = \"Scattering efficiency\",\n        title = \"BK7 glass spheres 4 micron diameter\",\n    )\nlines!(λ0 * 1000, qsca)\n\nfig","category":"page"},{"location":"#Water-Droplets","page":"Introduction","title":"Water Droplets","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"\"\"\"\nPlot the scattering cross section for 1 micron water droplets.\n\nThe plot shows the cross section as a function of wavelength.\n\"\"\"\n\nusing MieScattering\nusing CairoMakie\n\nnum = 100\nradius = 0.5 # in microns\nλ0 = LinRange(0.2, 1.2, 200)\nx = 2π * radius ./ λ0\n\n# from https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon - 24.0C\nm2 = @. 1.0 + 5.666959820E-1 / (1.0 - 5.084151894E-3 / λ0^2)\nm2 += @. 1.731900098E-1 / (1.0 - 1.818488474E-2 / λ0^2)\nm2 += @. 2.095951857E-2 / (1.0 - 2.625439472E-2 / λ0^2)\nm2 += @. 1.125228406E-1 / (1.0 - 1.073842352E1 / λ0^2)\nm = .√(m2)\nqext, qsca, qback, g = mie(m, x)\n\nfig = Figure()\nax =\n    fig[1, 1] = Axis(\n        fig,\n        xlabel = \"Wavelength [nm]\",\n        ylabel = \"Scattering Cross Section (µm²)\",\n        title = \"Water Droplets (1 µm diameter)\",\n    )\nlines!(λ0 * 1000, qsca)\n\nfig","category":"page"},{"location":"#Small-Gold-Spheres","page":"Introduction","title":"Small Gold Spheres","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"\"\"\"\nPlot the scattering cross section for 100nm gold spheres.\nThe resulting graph is as a function of wavelength.\n\"\"\"\n\nusing MieScattering\nusing CairoMakie\n\n# from https://refractiveindex.info/?shelf = main&book = Au&page = Johnson\n# wavelength in microns\nref_lam = [ 0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164,\n            0.2214, 0.2262, 0.2313, 0.2371, 0.2426, 0.2490, 0.2551, 0.2616,\n            0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107, 0.3204, 0.3315,\n            0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509,\n            0.4714, 0.4959, 0.5209, 0.5486, 0.5821, 0.6168, 0.6595, 0.7045,\n            0.7560, 0.8211, 0.8920, 0.9840, 1.0880, 1.2160, 1.3930, 1.6100, 1.9370]\n\nref_n = [   1.28, 1.32, 1.34, 1.33, 1.33, 1.30, 1.30, 1.30, 1.30, 1.31, 1.30,\n            1.32, 1.32, 1.33, 1.33, 1.35, 1.38, 1.43, 1.47, 1.49, 1.53, 1.53,\n            1.54, 1.48, 1.48, 1.50, 1.48, 1.46, 1.47, 1.46, 1.45, 1.38, 1.31,\n            1.04, 0.62, 0.43, 0.29, 0.21, 0.14, 0.13, 0.14, 0.16, 0.17, 0.22,\n            0.27, 0.35, 0.43, 0.56, 0.92]\n\nref_k = [   1.188, 1.203, 1.226, 1.251, 1.277, 1.304, 1.350, 1.387, 1.427,\n            1.460, 1.497, 1.536, 1.577, 1.631, 1.688, 1.749, 1.803, 1.847,\n            1.869, 1.878, 1.889, 1.893, 1.898, 1.883, 1.871, 1.866, 1.895,\n            1.933, 1.952, 1.958, 1.948, 1.914, 1.849, 1.833, 2.081, 2.455,\n            2.863, 3.272, 3.697, 4.103, 4.542, 5.083, 5.663, 6.350, 7.150,\n            8.145, 9.519, 11.21, 13.78]\n\nradius = 0.1 # in microns\nm = ref_n - im * ref_k\nx = 2π * radius ./ ref_lam\ncross_section_area = π * radius ^ 2\nmu_a = 4 * π * ref_k ./ ref_lam    # nm\nqext, qsca, qback, g = mie(m, x)\n\nsca_cross_section = qsca * cross_section_area\nabs_cross_section = (qext - qsca) * cross_section_area\n\nfig = Figure()\n\nax = fig[1,1] = Axis(fig, title=\"Gold Spheres 200nm diameter\")\nhidexdecorations!(ax)\nscatter!(ref_lam * 1000, ref_n, color=:blue, marker=:circle)\nscatter!(ref_lam * 1000, -ref_k, color=:red, marker=:rect)\ntext!(700, 1; text=\"real refractive index\", color=:blue)\ntext!(1100, -6; text=\"imaginary refractive index\", color=:red)\nylims!(-14,5)\n\nax = fig[2,1] = Axis(fig, ylabel=\"Absorption Depth [nm]\")\nscatter!(ref_lam * 1000, 1000 ./ mu_a, color=:blue, marker=:circle)\n\nax = fig[3,1] = Axis(fig, xlabel=\"Wavelength [nm]\", ylabel=\"Cross Section [µm²]\")\nscatter!(ref_lam * 1000, abs_cross_section, color=:blue, marker=:circle)\nscatter!(ref_lam * 1000, sca_cross_section, color=:red, marker=:rect)\ntext!(700, 0.01; text=\"absorption\", color=:blue)\ntext!(750, 0.1; text=\"scattering\", color=:red)\n\nfig","category":"page"}]
}
