module PulsarSearch
using StatsBase
include("stats.jl")

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
julia> using PulsarSearch

julia> z_n_binned(zeros(10), 2)
0.0

julia> z_n_binned(ones(10), 2) < 1e-30
true

julia> z_n_binned([10., 0., 0., 0., 0.], 2)
40.0
```
"""
function z_n_binned(profile::AbstractVector{}, n::Integer)::Float64
    total = sum(profile)
    N = length(profile)
    phase = range(0, stop = N - 1)  * (2pi/ N)

    if iszero(total)
        return 0.0
    end

    z = zero(Float64)
    for k in range(1, stop = n)
        s = 0.0
        c = 0.0
        for i in eachindex(profile)
            sk, ck = sincos(k * phase[i])
            s += profile[i] * sk
            c += profile[i] * ck
        end
        s = s^2
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
julia> using PulsarSearch

julia> z_n([10., 0., 0., 0., 0.], 2)
20.0
julia> z_n(ones(10), 2)
40.0
julia> z_n(Array([0.5]), 2)
0.0
```
"""
function z_n(phases::AbstractVector{T}, n::Integer)::T where {T<:AbstractFloat}
    N = length(phases)
    if N < 2
        return zero(T)
    end
    z = zero(T)
    twopiphase = 2 * pi * phases
    for k in range(1, stop = n)
        s = zero(T)
        c = zero(T)
        for ph in twopiphase
            sk, ck = sincos(k * ph)
            s += sk
            c += ck
        end
        s = s^2
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
    times::AbstractVector{T},
    n::Integer,
    fmin::Number,
    fmax::Number;
    oversample::Number = 2,
) where {T <: AbstractFloat}
    t0 = first(times)
    t1 = last(times)
    df = 1 / (t1 - t0) / oversample
    freqs = fmin:df:fmax
    N = length(freqs)
    stats = Vector{T}(undef, N)
    for i in eachindex(freqs)
        phases = times * freqs[i]
        phases .-= floor.(phases)
        stats[i] = z_n(phases, n)
    end
    return freqs, stats
end

"""
`z_n_search_hist(times, n, fmin, fmax [,oversample]) --> freqs, zsq_stat`

Calculate the ``Z^2_n`` statistics at trial frequencies in photon data.
Pre-bins the data using a histogram.
At the moment, it is _not_ faster. It will be after I add the 2-d hist
+ shift-and-add. But it allowed to test that `z_n_binned works` with
binned data.

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
function z_n_search_hist(
    times::AbstractVector{T},
    n::Integer,
    fmin::Number,
    fmax::Number;
    oversample::Number = 2,
    nbin::Integer = 16
) where {T <: AbstractFloat}
    t0 = first(times)
    t1 = last(times)
    df = 1 / (t1 - t0) / oversample
    freqs = fmin:df:fmax
    N = length(freqs)
    stats = Vector{T}(undef, N)
    for i in eachindex(freqs)
        phases = times * freqs[i]
        phases .-= floor.(phases)
        hist = fit(Histogram, phases, 0.0:(1/nbin):1.0)
        stats[i] = z_n_binned(hist.weights, n)
    end
    return freqs, stats
end


export z_n
export z_n_binned
export z_n_search
export z_n_search_hist
# Write your package code here.

end
