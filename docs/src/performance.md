# Performance in comparison to miepython

The following script shows how MieScattering.jl performs compared to miepython. The tests are performed by calling miepython using `PyCall`.
This might seem as an unfair comparison due to the (little) overhead of calling Python from Julia. However, if you are reading this you already want to use Julia for this task so resorting to PyCall might be your way to go anyway if there would be no native Julia solution.

The performance is compared vs the jit-ed and non-jit-ed version. Additionally the Julia version is used with the parameter `use_threads=false` to show the difference when using multi-threading in Julia.

```@example
using MieScattering
using BenchmarkTools
using PyCall
using Conda
Conda.pip_interop(true)
Conda.pip("install", "miepython")
using CairoMakie

miepython_jit = pyimport("miepython.miepython")
miepython = pyimport("miepython.miepython_nojit")

ntests = 6
m = 1.5
N = round.(Int, 10 .^ (LinRange(0, 3, ntests)))
N[1] = 2

result = zeros(ntests)
result_threads = zeros(ntests)
result_python = zeros(ntests)
result_python_jit = zeros(ntests)

for i = 1:ntests
    x = collect(LinRange(0.1, 20.0, N[i]))
    a = @benchmark mie($m, $x, use_threads = false)
    result[i] = median(a).time
    a = @benchmark mie($m, $x)
    result_threads[i] = median(a).time
    a = @benchmark miepython.mie($m, $x)
    result_python[i] = median(a).time
    a = @benchmark miepython_jit.mie($m, $x)
    result_python_jit[i] = median(a).time
end

imp = result_python_jit ./ result_threads

minimp, maximp = extrema(imp)

fig = Figure()
ax =
    fig[1, 1] = Axis(
        fig,
        yscale = log10,
        xscale = log10,
        xlabel = "Number of sphere sizes calculated",
        ylabel = "Execution Time [s]",
        title = "Julia improvement is $(round(minimp,digits=1))x to $(round(maximp,digits=1))x over numba (jit)",
    )
scatterlines!(N, result * 1e-9, label = "Julia (MieScattering.jl)")
scatterlines!(N, result_threads * 1e-9, label = "Julia Multi-Threading (MieScattering.jl)")
scatterlines!(N, result_python * 1e-9, label = "PyCall (miepython)")
scatterlines!(N, result_python_jit * 1e-9, label = "PyCall (miepython-jit)")
axislegend(ax; position = :lt)
fig

```