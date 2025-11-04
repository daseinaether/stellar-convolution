"""
Metallicity distribution function from COMPAS

DH0001_file
"""

import numpy as np
from scipy.stats import norm as NormDist


def compas_log_skew_normal_distribution_metallicity_distribution(
    redshifts,
    log_metallicity_centers,
    mu0,
    muz,
    sigma_0,
    sigma_z,
    alpha,
    global_logZ_distribution_min=-20,
    global_logZ_distribution_max=0,
    global_logZ_distribution_res=1000,
):
    """
    Calculate the distribution of metallicities at different redshifts using a log skew normal distribution
    the log-normal distribution is a special case of this log skew normal distribution distribution, and is retrieved by setting
    the skewness to zero (alpha = 0).

    Based on the method in Neijssel+19.

    Default values of mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0,
    retrieve the dP/dZ distribution used in Neijssel+19

    NOTE: This assumes that metallicities in COMPAS are drawn from a flat in log distribution!

    Args:
        max_redshift       --> [float]          max redshift for calculation
        redshift_step      --> [float]          step used in redshift calculation


        mu0    =  0.035    --> [float]           location (mean in normal) at redshift 0
        muz    = -0.25    --> [float]           redshift scaling/evolution of the location
        sigma_0 = 0.39     --> [float]          Scale (variance in normal) at redshift 0
        sigma_z = 0.00     --> [float]          redshift scaling of the scale (variance in normal)
        alpha   = 0.00    --> [float]          shape (skewness, alpha = 0 retrieves normal dist)

        min_logZ           --> [float]          Minimum logZ at which to calculate dPdlogZ (influences normalization)
        max_logZ           --> [float]          Maximum logZ at which to calculate dPdlogZ (influences normalization)
        step_logZ          --> [float]          Size of logZ steps to take in finding a Z range

    Returns:
        dPdlogZ            --> [2D float array] Probability of getting a particular logZ at a certain redshift
        p_draw_metallicity --> float            Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
    """

    #
    step_logZ = (
        global_logZ_distribution_max - global_logZ_distribution_min
    ) / global_logZ_distribution_res

    ##################################
    # create a range of metallicities (the x-values, or random variables)
    global_log_metallicities = np.arange(
        global_logZ_distribution_min,
        global_logZ_distribution_max + step_logZ,
        step_logZ,
    )

    ##################################
    # Log-Linear redshift dependence of sigma
    sigma = sigma_0 * 10 ** (sigma_z * redshifts)

    ##################################
    # Follow Langer & Norman 2007? in assuming that mean metallicities evolve in z as:
    mean_metallicities = mu0 * 10 ** (muz * redshifts)

    # Now we re-write the expected value of ou log-skew-normal to retrieve mu
    beta = alpha / (np.sqrt(1 + (alpha) ** 2))
    PHI = NormDist.cdf(beta * sigma)
    mu_metallicities = np.log(
        mean_metallicities / 2.0 * 1.0 / (np.exp(0.5 * sigma**2) * PHI)
    )

    ##################################
    # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
    dPdlogZ = (
        2.0
        / (sigma[:, np.newaxis])
        * NormDist.pdf(
            (global_log_metallicities - mu_metallicities[:, np.newaxis])
            / sigma[:, np.newaxis]
        )
        * NormDist.cdf(
            alpha
            * (global_log_metallicities - mu_metallicities[:, np.newaxis])
            / sigma[:, np.newaxis]
        )
    )

    ##################################
    # normalise the distribution over al metallicities
    norm = dPdlogZ.sum(axis=-1) * step_logZ
    dPdlogZ = dPdlogZ / norm[:, np.newaxis]

    ##################################
    # Select the metallicities that we have
    dPdLogZ_for_sampled_metallicities = dPdlogZ[
        :, np.digitize(log_metallicity_centers, global_log_metallicities) - 1
    ]

    # ##################################
    # # Calculate the dlogZ (stepsizes) values, adding one to the end.
    # dlogZ_sampled = np.diff(np.log(config["convolution_metallicity_bin_edges"]))

    # ##################################
    # # Calculate dP/dlogZ * dlogZ
    # dP = dPdLogZ_for_sampled_metallicities * dlogZ_sampled

    #
    return dPdLogZ_for_sampled_metallicities


def mean_metallicity(z, z0, alpha):
    """
    Function for mean metallicity
    """

    return z0 * np.power(10, alpha * z)


def mean_mu(z, z0, alpha, sigma):
    """
    Function to calculate mu
    """

    return np.log(mean_metallicity(z, z0, alpha)) - np.power(sigma, 2) / 2


def metallicity_distribution_lognormal(Z, z, z0, alpha, sigma):
    """
    Function to calculate the metallicity distribution on a grid of redshifts and metallicity according to a lognormal distribution.

    See Neijsel et al 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract) section 4.

    Function returns $dP(z)/dZ$
    """

    return (1 / (Z * sigma * np.power(2 * np.pi, 0.5))) * np.exp(
        -np.power(np.log(Z) - mean_mu(z, z0, alpha, sigma), 2)
        / (2 * np.power(sigma, 2))
    )


def metallicity_distribution_Neijsel19(log_metallicity_centers, redshifts):
    """
    Function to calculate the metallicity distribution fraction at a given redshift according to Neijsel et al 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)

    TODO: add kwargs so user can override anything
    """

    dPdlogZ = compas_log_skew_normal_distribution_metallicity_distribution(
        redshifts=redshifts,
        log_metallicity_centers=log_metallicity_centers,
        # metallicity distribution settings for Neijssel 2019
        mu0=0.035,
        muz=-0.23,
        sigma_0=0.39,
        sigma_z=0.0,
        alpha=0.0,
    )

    return dPdlogZ


def metallicity_distribution_vanSon2022(log_metallicity_centers, redshifts):
    """
    Function to calculate the metallicity distribution fraction as a function of redshift according to van Son et al. 2022

    TODO: add kwargs so user can override anything
    """

    dPdlogZ = compas_log_skew_normal_distribution_metallicity_distribution(
        redshifts=redshifts,
        log_metallicity_centers=log_metallicity_centers,
        # metallicity distribution settings for van Son 2021
        mu0=0.025,
        muz=-0.05,
        sigma_0=1.125,
        sigma_z=0.05,
        alpha=-1.77,
    )

    return dPdlogZ


def metallicity_distribution_dummy(constant=1):
    """
    Dummy metallicity distribution, based on the same functional form as neijsel et al 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)
    """

    # override redshift to have it not change
    return constant


if __name__ == "__main__":

    # compas_log_skew_normal_distribution_metallicity_distribution(redshifts=[0, 1, 2], logZmetallicity_centers=)

    dpdlogZ = metallicity_distribution_vanSon2022(
        log_metallicities=np.array([0.01, 0.005, 0.002, 0.001]),
        redshifts=np.array([0, 1, 2]),
    )
