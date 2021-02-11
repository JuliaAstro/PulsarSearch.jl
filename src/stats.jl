using Distributions


"""Equivalent gaussian sigma for small log-probability.

Return the equivalent gaussian sigma corresponding to the natural log of
the cumulative gaussian probability logp. In other words, return x, such
that Q(x) = p, where Q(x) is the cumulative normal distribution. This
version uses the rational approximation from Abramowitz and Stegun,
eqn 26.2.23, that claims to be precise to ~1e-4. Using the log(P) as input
gives a much extended range.

The parameters here are the result of a best-fit, with no physical meaning.

Translated from Scott Ransom's PRESTO
"""
function extended_equiv_gaussian_Nsigma(logp::Number)
    t = sqrt(-2.0 * logp)
    num = 2.515517 + t * (0.802853 + t * 0.010328)
    denom = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308))
    return t - num / denom
end


"""Number of Gaussian sigmas corresponding to tail log-probability.

This function computes the value of the characteristic function of a
standard Gaussian distribution for the tail probability equivalent to the
provided p-value, and turns this value into units of standard deviations
away from the Gaussian mean. This allows the user to make a statement
about the signal such as “I detected this pulsation at 4.1 sigma

The example values below are obtained by brute-force integrating the
Gaussian probability density function using the mpmath library
between Nsigma and +inf.

# Examples
```jldoctest
julia> using psrsearch

julia> pvalues = Array([0.15865525393145707, 0.0013498980316301035, 9.865877e-10, 6.22096e-16, 3.0567e-138]);

julia> log_pvalues = log.(pvalues);

julia> sigmas = Array([1, 3, 6, 8, 25]);

julia> all(isapprox(equivalent_gaussian_Nsigma_from_logp.(log_pvalues), sigmas; atol=0.1))
true
```
"""
function equivalent_gaussian_Nsigma_from_logp(logp)

    if logp < -300
        # print("Extended")
        return extended_equiv_gaussian_Nsigma(logp)
    end
    d = Normal()
    return cquantile(d, exp(logp))
end

export equivalent_gaussian_Nsigma_from_logp

"""Number of Gaussian sigmas corresponding to tail probability.

This function computes the value of the characteristic function of a
standard Gaussian distribution for the tail probability equivalent to the
provided p-value, and turns this value into units of standard deviations
away from the Gaussian mean. This allows the user to make a statement
about the signal such as “I detected this pulsation at 4.1 sigma
"""
function equivalent_gaussian_Nsigma(p)
    return equivalent_gaussian_Nsigma_from_logp(log(p))
end

export equivalent_gaussian_Nsigma

"""Asymptotic natural log of incomplete gamma function.

Return the natural log of the incomplete gamma function in
its asymptotic limit as z->infty.  This is from Abramowitz
and Stegun eqn 6.5.32.

Translated from Scott Ransom's PRESTO



# Examples

```jldoctest
julia> using psrsearch

julia> pvalues = Array([0.15865525393145707, 0.0013498980316301035]);

julia> sigmas = Array([1, 3]);

julia> all(isapprox(equivalent_gaussian_Nsigma.(pvalues), sigmas; atol=0.1))
true

```

"""
function log_asymptotic_incomplete_gamma(a, z)
    x = 1.0
    newxpart = 1.0
    term = 1.0
    ii = 1

    while (abs(newxpart) > 1e-15)
        term *= (a - ii)
        newxpart = term / z^ii
        x += newxpart
        ii += 1
    end

    return (a - 1.0) * log(z) - z + log(x)
end

export log_asymptotic_incomplete_gamma

"""Natural log of the Gamma function in its asymptotic limit.

Return the natural log of the gamma function in its asymptotic limit
as z->infty.  This is from Abramowitz and Stegun eqn 6.1.41.

Translated from Scott Ransom's PRESTO
"""
function log_asymptotic_gamma(z)
    half_log_twopi = 0.91893853320467267  # (1/2)*log(2*pi)
    one_twelfth = 8.3333333333333333333333e-2
    one_degree = 2.7777777777777777777778e-3  # 1 / 360
    one_over_1680 = 5.9523809523809529e-4
    one_over_1260 = 7.9365079365079365079365e-4
    x = (z - 0.5) * log(z) - z + half_log_twopi
    y = 1.0 / (z * z)
    x +=
        (
            ((-one_over_1680 * y + one_over_1260) * y - one_degree) * y +
            one_twelfth
        ) / z
    return x
end

"""Log survival function of the chi-squared distribution.

# Examples

```jldoctest chi2_logp
julia> using psrsearch;

julia> using Distributions;

julia> chi2 = 31;

julia> d = Chisq(2);

julia> isapprox(chi2_logp(chi2, 2), logccdf(d, chi2), atol=0.1)
true

julia> chi2 = Array([5, 32]);

julia> all(isapprox.(chi2_logp.(chi2, 2), logccdf.(d, chi2), atol=0.1))
true
```
"""
function chi2_logp(chi2, dof)

    # If very large reduced chi squared, use approximation. This is an
    # eyeballed limit parameter space where the difference between the
    # approximation and the scipy version is tiny, but above which the scipy
    # version starts failing.
    if (chi2 / dof > 15.0) || ((dof > 150) && (chi2 / dof > 6.0))
        return log_asymptotic_incomplete_gamma(0.5 * dof, 0.5 * chi2) -
               log_asymptotic_gamma(0.5 * dof)
    end
    d = Chisq(dof)
    return logccdf(d, chi2)
end

export chi2_logp

"""Calculate a multi-trial p-value from the log of a single-trial one.

This allows to work around Numba's limitation on longdoubles, a way to
vectorize the computation when we need longdouble precision.

Parameters
----------
logp1 : float
    The natural logarithm of the significance at which we reject the null
    hypothesis on each single trial.
n : int
    The number of trials

Returns
-------
logpn : float
    The log of the significance at which we reject the null hypothesis
    after multiple trials
"""
function logp_multitrial_from_single_logp(logp1, n)
    # If the the probability is very small (p1 * n) < 1e-6, use Bonferroni
    # approximation.
    logn = log(n)
    if logp1 + logn < -7
        return logp1 + logn
    end

    return log(1 - (1 - exp(logp1))^n)
end

export logp_multitrial_from_single_logp

raw"""Calculate a multi-trial p-value from a single-trial one.

Calling *p* the probability of a single success, the Binomial
distributions says that the probability *at least* one outcome
in n trials is

``P(k\geq 1) = \sum_{k\geq 1} \binom{n}{k} p^k (1-p)^{(n-k)}``

or more simply, using P(k ≥ 0) = 1

P(k\geq 1) = 1 - \binom{n}{0} (1-p)^n = 1 - (1-p)^n``


Parameters
----------
p1 : float
    The significance at which we reject the null hypothesis on
    each single trial.
n : int
    The number of trials

Returns
-------
pn : float
    The significance at which we reject the null hypothesis
    after multiple trials
"""
function p_multitrial_from_single_trial(p1, n)
    logpn = logp_multitrial_from_single_logp(log(p1), n)
    return exp(logpn)
end

export p_multitrial_from_single_trial

"""Calculate a multi-trial p-value from the log of a single-trial one.

This allows to work around Numba's limitation on longdoubles, a way to
vectorize the computation when we need longdouble precision.

Parameters
----------
logpn : float
    The natural logarithm of the significance at which we want to reject
    the null hypothesis after multiple trials
n : int
    The number of trials

Returns
-------
logp1 : float
    The log of the significance at which we reject the null hypothesis on
    each single trial.
"""
function logp_single_trial_from_logp_multitrial(logpn, n)
    logn = log(n)
    # If the the probability is very small, use Bonferroni approximation.
    if logpn < -7
        return logpn - logn
    end
    # Numerical errors arise when pn is very close to 1.
    p1 = 1 - (1 - exp(logpn))^(1 / n)
    return log(p1)
end

export logp_single_trial_from_logp_multitrial

raw"""Calculate the single-trial p-value from a total p-value

Let us say that we want to reject a null hypothesis at the
``pn`` level, after executing ``n`` different measurements.
This might be the case because, e.g., we
want to have a 1% probability of detecting a signal in an
entire power spectrum, and we need to correct the detection
level accordingly.

The typical procedure is dividing the initial probability
(often called _epsilon_) by the number of trials. This is
called the Bonferroni correction and it is often a good
approximation, when ``pn`` is low: ``p1 = pn / n``.

However, if ``pn`` is close to 1, this approximation gives
incorrect results.

Here we calculate this probability by inverting the Binomial
problem. Given that (see ``p_multitrial_from_single_trial``)
the probability of getting more than one hit in n trials,
given the single-trial probability *p*, is

``P (k \geq 1) =  1 - (1 - p)^n``

we get the single trial probability from the multi-trial one
from

``p = 1 - (1 - P)^{(1/n)}``

This is also known as Šidák correction.

Parameters
----------
pn : float
    The significance at which we want to reject the null
    hypothesis after multiple trials
n : int
    The number of trials

Returns
-------
p1 : float
    The significance at which we reject the null hypothesis on
    each single trial.
"""
function p_single_trial_from_p_multitrial(pn, n)
    logp = logp_single_trial_from_logp_multitrial(log(pn), n)

    return exp(logp)
end

export p_single_trial_from_p_multitrial

"""Calculate the probability of a certain folded profile, due to noise.

Parameters
----------
z2 : float
    A ``Z^2_n`` statistics value
n : int, default 2
    The ``n`` in ``Z^2_n`` (number of harmonics, including the fundamental)

Other Parameters
----------------
ntrial : int
    The number of trials executed to find this profile
n_summed_spectra : int
    Number of ``Z_2^n`` periodograms that were averaged to obtain z2

Returns
-------
p : float
    The probability that the ``Z^2_n`` value has been produced by noise
"""
function z2_n_probability(z2, n; ntrial = 1, n_summed_spectra = 1)
    d = Chisq(2 * n * n_summed_spectra)
    epsilon_1 = ccdf(d, z2 * n_summed_spectra)
    epsilon = p_multitrial_from_single_trial(epsilon_1, ntrial)
    return epsilon
end

export z2_n_probability

"""Calculate the probability of a certain folded profile, due to noise.

Parameters
----------
z2 : float
    A ``Z^2_n`` statistics value
n : int, default 2
    The ``n`` in ``Z^2_n`` (number of harmonics, including the fundamental)

Other Parameters
----------------
ntrial : int
    The number of trials executed to find this profile
n_summed_spectra : int
    Number of ``Z_2^n`` periodograms that were averaged to obtain z2

Returns
-------
p : float
    The probability that the ``Z^2_n`` value has been produced by noise
"""
function z2_n_logprobability(z2, n; ntrial = 1, n_summed_spectra = 1)

    epsilon_1 = chi2_logp(z2 * n_summed_spectra, 2 * n * n_summed_spectra)
    epsilon = logp_multitrial_from_single_logp(epsilon_1, ntrial)
    return epsilon
end

export z2_n_logprobability

"""Return the detection level for the ``Z^2_n`` statistics.

See Buccheri et al. (1983), Bendat and Piersol (1971).

Parameters
----------
n : int, default 2
    The ``n`` in ``Z^2_n`` (number of harmonics, including the fundamental)
epsilon : float, default 0.01
    The fractional probability that the signal has been produced by noise

Other Parameters
----------------
ntrial : int
    The number of trials executed to find this profile
n_summed_spectra : int
    Number of Z_2^n periodograms that are being averaged

Returns
-------
detlev : float
    The epoch folding statistics corresponding to a probability
    epsilon * 100 % that the signal has been produced by noise
"""
function z2_n_detection_level(
    n = 2,
    epsilon = 0.01;
    ntrial = 1,
    n_summed_spectra = 1,
)
    epsilon = p_single_trial_from_p_multitrial(epsilon, ntrial)
    d = Chisq(2 * n_summed_spectra * n)
    retlev = cquantile(d, epsilon) / (n_summed_spectra)
    return retlev
end

export z2_n_detection_level
