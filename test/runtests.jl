using MieScattering
using Test

@testset "MieScattering.jl" begin
    include("suites/nonabsorbing.jl")
    include("suites/absorbing.jl")
    include("suites/small.jl")
    include("suites/perfectlyreflecting.jl")
    include("suites/anglescattering.jl")
    include("suites/phasematrix.jl")
    include("suites/notebook.jl")
end
