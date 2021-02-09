using psrsearch
using Test
using Distributions

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
