__precompile__()

module psrsearch
using Distributions


function z_n_binned(profile::Array{T,1}, n::Integer)::T where {T<:AbstractFloat}
    total = sum(profile)
    N = length(profile)
    phase = Array(range(0, stop = N - 1) / N) * 2 * pi

    if total == 0
        return 0
    end

    z = 0
    for k in range(1, stop = n)
        s = 0
        for i in range(1, stop = N)
            s += profile[i] * sin(k * phase[i])
        end
        s = s^2
        c = 0
        for i in range(1, stop = N)
            c += profile[i] * cos(k * phase[i])
        end
        c = c^2
        z += c + s
    end
    return z * 2 / total
end


function z_n(phases::Array{T,1}, n::Integer)::T where {T<:AbstractFloat}
    N = length(phases)
    z = 0
    twopiphase = 2 * pi * phases
    for k in range(1, stop = n)
        s = 0
        for ph in twopiphase
            s += sin(k * ph)
        end
        s = s^2
        c = 0
        for ph in twopiphase
            c += cos(k * ph)
        end
        c = c^2
        z += c + s
    end
    return z * 2 / N
end


function z_n_search(
    times::Array{T,1},
    n::Integer,
    fmin::Number,
    fmax::Number;
    oversample::Number = 2,
) where {T <: AbstractFloat}
    t0 = first(times)
    t1 = last(times)
    df = 1 / (t1 - t0) / oversample
    freqs = [fmin:df:fmax;]
    N = length(freqs)
    stats = zeros(N)
    for i in range(1, stop = N)
        phases = times * freqs[i]
        phases -= floor.(phases)
        stats[i] = z_n(phases, n)
    end
    return freqs, stats
end

export z_n
export z_n_binned
export z_n_search
# Write your package code here.

end
