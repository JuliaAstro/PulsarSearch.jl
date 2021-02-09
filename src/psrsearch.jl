__precompile__()

module psrsearch
using Distributions

"""
`z_n_binned(profile, n) --> zsq`

``Z^2_n`` statistic for pulse profiles from binned events

See Bachetti+2021, arXiv:2012.11397

Parameters
----------
`profile` : array of floats
    The folded pulse profile (containing the number of
    photons falling in each pulse bin)

`n` : int
    Number of harmonics, including the fundamental

Returns
-------
`zsq` : float
    The value of the statistic

# Examples
```jldoctest
julia> using psrsearch

julia> z_n_binned(zeros(10), 2)
0.0

julia> z_n_binned(ones(10), 2) < 1e-30
true

julia> z_n_binned([10., 0., 0., 0., 0.], 2)
40.0
```
"""
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

"""
`z_n(phases, n) --> zsq`

``Z^2_n`` statistics, à la Buccheri+83, A&A, 128, 245, eq. 2.

Parameters
----------
`phase` : array of floats
    The phases of the events

`n` : int, default 2
    Number of harmonics, including the fundamental

Returns
-------
`zsq` : float
    The ``Z^2_n`` statistic

# Examples
```jldoctest
julia> using psrsearch

julia> z_n([10., 0., 0., 0., 0.], 2)
20.0
julia> z_n(ones(10), 2)
40.0
julia> z_n(Array([0.5]), 2)
0.0
```
"""
function z_n(phases::Array{T,1}, n::Integer)::T where {T<:AbstractFloat}
    N = length(phases)
    if N < 2
        return 0
    end
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

"""
`z_n_search(times, n, fmin, fmax [,oversample]) --> freqs, zsq_stat`

Calculate the ``Z^2_n`` statistics at trial frequencies in photon data.

Parameters
----------
`times` : array-like
    the event arrival times

`n` : int
    the number of harmonics in ``Z^2_n``

Other Parameters
----------------
`fmin` : float
    minimum pulse frequency to search

`fmax` : float
    maximum pulse frequency to search

`oversample` : float
    Oversampling factor with respect to the usual 1/T/n rule

Returns
-------
`fgrid` : array-like
    frequency grid of the epoch folding periodogram

`zsq_stat` : array-like
    the Z^2_n statistics corresponding to each frequency bin.
"""
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
