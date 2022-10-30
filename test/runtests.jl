using MieScattering
using Test

@testset "MieScattering.jl" begin
    @testset "Non-Absorbing" begin
        m = 1.55
        λ0 = 0.6328
        radius = 0.525
        x = 2π*radius / λ0
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 3.10543 atol=0.00001
        @test qsca ≈ 3.10543 atol=0.00001
        @test qback ≈ 2.92534 atol=0.00001
        @test g ≈ 0.63314 atol=0.00001


        # MIEV0 Test Case 5
        m = 0.75 + 0.0im
        x = 0.099
        s1 = 1.81756e-8 - 1.64810e-4*im
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000007 atol=1e-6
        @test g ≈ 0.001448 atol=1e-6
        @test qback ≈ G atol=1e-6

        # MIEV0 Test Case 6
        m = 0.75 + 0.0im
        x = 0.101
        s1 = 2.04875e-8 - im*1.74965e-4
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000008 atol=1e-6
        @test g ≈ 0.001507 atol=1e-6
        @test qback ≈ G atol=1e-6

        # MIEV0 Test Case 7
        m = 0.75 + 0.0im
        x = 10.0
        s1 = -1.07857 - im*3.60881e-2
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.232265 atol=1e-6
        @test g ≈ 0.896473 atol=1e-6
        @test qback ≈ G atol=1e-6

        # MIEV0 Test Case 8
        m = 0.75 + 0.0im
        x = 1000.0
        s1 = 1.70578e1 - im*4.84251e2
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 1.997908 atol=1e-6
        @test g ≈ 0.844944 atol=1e-6
        @test qback ≈ G atol=1e-6

        # OLD MIEV0 Test Case 1
        m = 1.5 + 0.0im
        x = 10
        s1 = 4.322 - im*4.868
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.882 atol=1e-4
        @test qback ≈ G atol=1e-4

        # OLD MIEV0 Test Case 2
        m = 1.5 + 0.0im
        x = 100
        s1 = 4.077e1 + im*5.175e1
        G = abs2(2*s1/x)
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.0944 atol=1e-4
        @test qback ≈ G atol=1e-4

        # OLD MIEV0 Test Case 3
        m = 1.5 + 0.0im
        x = 1000
        G = 4*2.576e6/x^2
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.0139 atol=1e-4
        @test qback ≈ G atol=1e-3

        # OLD MIEV0 Test Case 4
        m = 1.5 + 0.0im
        x = 5000
        G = 4*2.378e8/x^2
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.0086 atol=1e-4
        @test qback ≈ G atol=3e-3

        m = 1.55 - im*0.1
        λ0 = 0.6328
        radius = 0.525
        x = 2π*radius/λ0
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 2.86165188243 atol=1e-7
        @test qsca ≈ 1.66424911991 atol=1e-7
        @test qback ≈ 0.20599534080 atol=1e-7
        @test g ≈ 0.80128972639 atol=1e-7
    end

    @testset "Absorbing" begin
        # MIEV0 Test Case 9
        m = 1.33 - im*0.00001
        x = 1.0
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.093923 atol=1e-6
        @test g ≈ 0.184517 atol=1e-6

        # MIEV0 Test Case 10
        m = 1.33 - im*0.00001
        x = 100.0
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.096594 atol=1e-6
        @test g ≈ 0.868959 atol=1e-6

        # MIEV0 Test Case 11
        m = 1.33 - im*0.00001
        x = 10000.0
        qext, qsca, qback, g = mie(m, x)
        @test g ≈ 0.907840 atol=1e-6
        @test qsca ≈ 1.723857 atol=1e-6

        # MIEV0 Test Case 12
        m = 1.5 - im
        x = 0.055
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000011 atol=1e-6
        @test g ≈ 0.000491 atol=1e-6

        # MIEV0 Test Case 13
        m = 1.5 - im
        x = 0.056
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000012 atol=1e-6
        @test g ≈ 0.000509 atol=1e-6

        # MIEV0 Test Case 14
        m = 1.5 - im
        x = 1
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.6634538 atol=1e-6
        @test g ≈ 0.192136 atol=1e-6

        # MIEV0 Test Case 15
        m = 1.5 - im
        x = 100.0
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 1.283697 atol=1e-3
        @test qext ≈ 2.097502 atol=1e-2
        @test g ≈ 0.850252 atol=1e-3

        # MIEV0 Test Case 16
        m = 1.5 - im
        x = 10000
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 1.236575 atol=1e-6
        @test qext ≈ 2.004368 atol=1e-6
        @test g ≈ 0.846309 atol=1e-6

        # MIEV0 Test Case 17
        m = 10.0 - im*10.0
        x = 1
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.049405 atol=1e-6
        @test g ≈ -0.110664 atol=1e-6

        # MIEV0 Test Case 18
        m = 10.0 - im*10.0
        x = 100
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 1.836785 atol=1e-6
        @test g ≈ 0.556215 atol=1e-6

        # MIEV0 Test Case 19
        m = 10.0 - im*10.0
        x = 10000
        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 1.795393 atol=1e-6
        @test g ≈ 0.548194 atol=1e-6

        m = 1.5 - im*0.5
        x = 2.5
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 2.562873497454734 atol=1e-7
        @test qsca ≈ 1.097071819088392 atol=1e-7
        @test qback ≈ 0.123586468179818 atol=1e-7
        @test g ≈ 0.748905978948507 atol=1e-7
    end

    @testset "Perfectly reflecting" begin
        # MIEV0 Test Case 0
        m = 0
        x = 0.001

        qext, qsca, qback, g = mie(m, x)
        @test isapprox(qsca, 3.3333e-12, atol=1e-13)

        # MIEV0 Test Case 1
        m = 0
        x = 0.099

        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000321 atol=1e-4
        @test g ≈ -0.397357 atol=1e-3

        # MIEV0 Test Case 2
        m = 0
        x = 0.101

        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 0.000348 atol=1e-6
        @test g ≈ -0.397262 atol=1e-6

        # MIEV0 Test Case 3
        m = 0
        x = 100

        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.008102 atol=1e-6
        @test g ≈ 0.500926 atol=1e-6

        # MIEV0 Test Case 4
        m = 0
        x = 10000

        qext, qsca, qback, g = mie(m, x)
        @test qsca ≈ 2.000289 atol=1e-6
        @test g ≈ 0.500070 atol=1e-6
    end

   @testset "Small" begin
        # MIEV0 Test Case 5
        m = 0.75
        x = 0.099

        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000007 atol=1e-6
        @test g ≈ 0.001448 atol=1e-6

        # MIEV0 Test Case 6
        m = 0.75
        x = 0.101

        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000008 atol=1e-6
        @test g ≈ 0.001507 atol=1e-6

        m = 1.5 - im
        x = 0.055
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.101491 atol=1e-6
        @test g ≈ 0.000491 atol=1e-6
        x = 0.056
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.103347 atol=1e-6
        @test g ≈ 0.000509 atol=1e-6

        m = 1e-10 - im*1e10
        x = 0.099
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000321 atol=1e-6
        @test g ≈ -0.397357 atol=1e-4
        x = 0.101
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000348 atol=1e-6
        @test g ≈ -0.397262 atol=1e-6

        m = 0 - im*1e10
        x = 0.099
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000321 atol=1e-6
        @test g ≈ -0.397357 atol=1e-4
        x = 0.101
        qext, qsca, qback, g = mie(m, x)
        @test qext ≈ 0.000348 atol=1e-6
        @test g ≈ -0.397262 atol=1e-4
   end

   @testset "AngleScattering" begin
        x = 1.0
        m = 1.5 - im*1.0
        θ = 0:30:180
        μ = cosd.(θ)

        qext, qsca, qback, g = mie(m, x)
        S1, S2 = mie_S1_S2(m, x, μ)
        S1 *= √(π*x^2*qext)
        S2 *= √(π*x^2*qext)

        @test real(S1[1]) ≈ 0.584080 atol=1e-6
        @test imag(S1[1]) ≈ 0.190515 atol=1e-6
        @test real(S2[1]) ≈ 0.584080 atol=1e-6
        @test imag(S2[1]) ≈ 0.190515 atol=1e-6

        @test real(S1[2]) ≈ 0.565702 atol=1e-6
        @test imag(S1[2]) ≈ 0.187200 atol=1e-6
        @test real(S2[2]) ≈ 0.500161 atol=1e-6
        @test imag(S2[2]) ≈ 0.145611 atol=1e-6

        @test real(S1[3]) ≈ 0.517525 atol=1e-6
        @test imag(S1[3]) ≈ 0.178443 atol=1e-6
        @test real(S2[3]) ≈ 0.287964 atol=1e-6
        @test imag(S2[3]) ≈ 0.041054 atol=1e-6

        @test real(S1[4]) ≈ 0.456340 atol=1e-6
        @test imag(S1[4]) ≈ 0.167167 atol=1e-6
        @test real(S2[4]) ≈ 0.0362285 atol=1e-6
        @test imag(S2[4]) ≈ -0.0618265 atol=1e-6

        @test real(S1[5]) ≈ 0.400212 atol=1e-6
        @test imag(S1[5]) ≈ 0.156643 atol=1e-6
        @test real(S2[5]) ≈ -0.174875 atol=1e-6
        @test imag(S2[5]) ≈ -0.122959 atol=1e-6

        @test real(S1[6]) ≈ 0.362157 atol=1e-6
        @test imag(S1[6]) ≈ 0.149391 atol=1e-6
        @test real(S2[6]) ≈ -0.305682 atol=1e-6
        @test imag(S2[6]) ≈ -0.143846 atol=1e-6

        @test real(S1[7]) ≈ 0.348844 atol=1e-6
        @test imag(S1[7]) ≈ 0.146829 atol=1e-6
        @test real(S2[7]) ≈ -0.348844 atol=1e-6
        @test imag(S2[7]) ≈ -0.146829 atol=1e-6
   end

   @testset "Notebook" begin
        N = 500
        m = 1.5
        x = LinRange(0.1, 20, N)
        qext, qsca, qback, g = mie(m, x)

        @test qsca[1] ≈ 2.3084093592198083e-05 atol=1e-6
        @test qsca[101] ≈ 4.105960809066763 atol=1e-6
        @test qsca[201] ≈ 1.9947867190110644 atol=1e-6
        @test qsca[301] ≈ 2.4652591512196405 atol=1e-6
        @test qsca[401] ≈ 2.472171798724846 atol=1e-6
        @test qsca[500] ≈ 2.03583698038088 atol=1e-6

        N = 500
        m = 1.5
        ρ = LinRange(0.1, 20, N)
        x15 = ρ/2/(m-1)
        qext, scal5, qback, g = mie(m, x15)

        m = 1.1
        x11 = ρ/2/(m-1)
        qext, scal1, qback, g = mie(m, x11)

        @test scal1[1] ≈ 0.0006616369953521216 atol=1e-6
        @test scal1[100] ≈ 3.449616595439377 atol=1e-6
        @test scal1[200] ≈ 1.6837703285684387 atol=1e-6
        @test scal1[300] ≈ 2.3167184401740495 atol=1e-6
        @test scal1[400] ≈ 2.218210809017406 atol=1e-6
        @test scal1[500] ≈ 1.876467571615533 atol=1e-6

        @test scal5[1] ≈ 2.3084093592198083e-05 atol=1e-6
        @test scal5[100] ≈ 4.07295075914037 atol=1e-6
        @test scal5[200] ≈ 1.8857586341949146 atol=1e-6
        @test scal5[300] ≈ 2.464763930426085 atol=1e-6
        @test scal5[400] ≈ 2.430569030744473 atol=1e-6
        @test scal5[500] ≈ 2.03583698038088 atol=1e-6
   end
end
