using MieScattering
using BenchmarkTools

ntests = 6
m=1.5
N = round.(Int, 10 .^(LinRange(0, 3, ntests)))
N[1] = 2

result = zeros(ntests)

for i=1:ntests
    x = LinRange(0.1, 20.0, N[i])
    a = @benchmark mie($m, $x)
    result[i] = median(a).time
end
