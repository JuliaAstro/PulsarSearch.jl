using psrsearch
using Test
using Distributions
using Documenter

doctest(psrsearch)

@testset "psrsearch.jl" begin
    @testset "Basic" begin
        flatarr = Array([0.5])
        @test z_n(flatarr, 2) == 0
        flatarr = Array(range(1, stop = 10) / 10)
        @test z_n(flatarr, 2) < 1e-10
        @test z_n(ones(10), 2) == 40
    end

    @testset "Binned" begin
        flatarr = zeros(10)
        @test z_n_binned(flatarr, 2) == 0
        flatarr = ones(10)
        @test z_n_binned(flatarr, 2) < 1e-10
        arr = zeros(10)
        arr[1] = 10
        @test z_n(arr, 2) == 40
    end

    @testset "Search" begin
        phases = rand(Normal(0.5, 0.1), 1000)
        pulse_no = rand(Uniform(0, 1000), 1000)
        pulse_no = floor.(pulse_no)
        f = 1.123
        times = sort((phases + pulse_no) / f)
        freqs, stats = z_n_search(times, 2, 1.0, 1.5, oversample = 4.0)
        maxind = argmax(stats)
        @test abs(freqs[maxind] - f) < 1e-3
    end
end

@testset "Stat" begin
    @testset "single_from_multi $ntrial" for ntrial in [1, 10, 100, 1000, 10000, 100000]
        epsilon_1 = 0.00000001
        epsilon_n = p_multitrial_from_single_trial(epsilon_1, ntrial)
        epsilon_1_corr = p_single_trial_from_p_multitrial(epsilon_n, ntrial)

        @test isapprox(epsilon_1_corr, epsilon_1; rtol=1e-2)
    end
    @testset "Zn Det Lev" begin
        @test isapprox(z2_n_detection_level(2), 13.276704135987625)
        epsilon_corr = p_single_trial_from_p_multitrial(0.01, 2)
        @test isapprox(z2_n_detection_level(4, 0.01, ntrial=2),
                        z2_n_detection_level(4, epsilon_corr))

    end
    @testset "Zn Det Lev ntrial $ntrial" for ntrial in [1, 10, 100, 1000, 100000]
        detlev = z2_n_detection_level(2, 0.1, ntrial=ntrial)
        @test isapprox(z2_n_probability(detlev, 2, ntrial=ntrial), 0.1)
    end
end
