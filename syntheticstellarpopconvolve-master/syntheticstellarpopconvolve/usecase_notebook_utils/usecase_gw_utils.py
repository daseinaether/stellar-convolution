import os

import deepdish as dd
import matplotlib.patches as patches
import numpy as np

KAPPA_DEFAULT = 2.9
REDSHIFT_DEFAULT = 0.2


#######################
# primary mass
def get_data_powerlaw_peak_primary_mass(data_root, redshift, limits, kappa):
    """
    Routine to get the data for the powerlaw peak estimates
    """

    # Get file and set limits
    mass_PP_path = os.path.join(
        data_root, "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_mass_data.h5"
    )

    # Create mass grid
    mass_1 = np.linspace(2, 100, 1000)
    mass_ratio = np.linspace(0.1, 1, 500)

    # load in the traces.
    # Each entry in lines is p(m1 | Lambda_i) or p(q | Lambda_i)
    # where Lambda_i is a single draw from the hyperposterior
    # The ppd is a 2D object defined in m1 and q
    with open(mass_PP_path, "r") as _data:
        _data = dd.io.load(mass_PP_path)
        lines = _data["lines"]
        ppd = _data["ppd"]

    # Set redshift scaling multiplication factor
    redshift_scaling_factor = (1 + redshift) ** (kappa)

    # marginalize over q to get the ppd in terms of m1 only
    mass_1_ppd = np.trapz(ppd, mass_ratio, axis=0) * redshift_scaling_factor
    CI_down = (
        np.percentile(lines["mass_1"], limits[0], axis=0) * redshift_scaling_factor
    )
    CI_up = np.percentile(lines["mass_1"], limits[1], axis=0) * redshift_scaling_factor

    return_dict = {
        "mass_1": mass_1,
        "mass_1_ppd": mass_1_ppd,
        "mass_1_lines": lines["mass_1"],
        "CI_up": CI_up,
        "CI_down": CI_down,
    }

    return return_dict


def add_primary_mass_distribution_to_figure(
    fig, ax, mass_1, mass_1_ppd, CI_down, CI_up, label=None, fill_between_kwargs=None
):
    """
    Function to add the distribution to the figure
    """

    if fill_between_kwargs is None:
        fill_between_kwargs = {}

    # plot the PPD as a solid line
    ax.semilogy(
        mass_1,
        mass_1_ppd,
        label=label,
        alpha=0.75,
        **fill_between_kwargs,
    )

    # plot the CIs as a filled interval
    ax.fill_between(
        mass_1,
        CI_down,
        CI_up,
        alpha=0.5,
        # label=label,
        **fill_between_kwargs,
    )

    return fig, ax


def add_text_to_primary_mass_distribution_plot(fig, ax, redshift, kappa):
    """
    Function to add text to the primary mass distribution plot
    """

    # Set redshift scaling multiplication factor
    redshift_scaling_factor = (1 + redshift) ** (kappa)

    # with plt.xkcd():
    ax.text(
        20,
        0.5 * 10**-1 * redshift_scaling_factor,
        "GWTC-3\nPL+Peak model",
        # "GWTC-3\nPL+Peak model\n(z={} \kappa={})".format(redshift, kappa),
        rotation=-30,
        horizontalalignment="right",
        verticalalignment="center",
        fontsize=20,
    )

    #
    a = patches.FancyArrowPatch(
        (20, 0.5 * 10**-1 * redshift_scaling_factor),
        (22, 0.65 * 10**-1 * redshift_scaling_factor),
        connectionstyle="arc3,rad=.4",
        arrowstyle=patches.ArrowStyle.Fancy(head_length=8, head_width=12, tail_width=2),
        color="k",
    )
    ax.add_patch(a)

    return fig, ax


def add_confidence_interval_powerlaw_peak_primary_mass(
    fig,
    ax,
    data_root,
    add_text=False,
    label=None,
    fill_between_kwargs=None,
    limits=[5, 95],
    redshift=REDSHIFT_DEFAULT,
    kappa=KAPPA_DEFAULT,
):
    """
    Function to add the confidence interval for the powerlaw + peak model of the GWTC03b data release for the primary mass distribution

    data_root has to contain the file "o1o2o3_mass_c_iid_mag_iid_tilt_powerlaw_redshift_mass_data.h5"
    """

    #
    if fill_between_kwargs is None:
        fill_between_kwargs = {}

    # get data
    data_powerlaw_peak_primary_mass = get_data_powerlaw_peak_primary_mass(
        data_root=data_root, redshift=redshift, limits=limits, kappa=kappa
    )

    # unpack
    mass_1 = data_powerlaw_peak_primary_mass["mass_1"]
    mass_1_ppd = data_powerlaw_peak_primary_mass["mass_1_ppd"]
    CI_up = data_powerlaw_peak_primary_mass["CI_up"]
    CI_down = data_powerlaw_peak_primary_mass["CI_down"]

    fig, ax = add_primary_mass_distribution_to_figure(
        fig=fig,
        ax=ax,
        mass_1=mass_1,
        mass_1_ppd=mass_1_ppd,
        CI_down=CI_down,
        CI_up=CI_up,
        label=label,
        fill_between_kwargs=fill_between_kwargs,
    )

    # Add text to plot
    if add_text:
        fig, ax = add_text_to_primary_mass_distribution_plot(
            fig=fig, ax=ax, redshift=redshift, kappa=kappa
        )

    return fig, ax


#########################


def get_median_percentiles(value_array):
    """
    Function to get the median and the percentiles from the data
    """

    result_dict = {}

    result_dict["median"] = np.percentile(value_array, [50], axis=0)

    result_dict["90%_CI"] = np.percentile(value_array, [5, 95], axis=0)

    # result_dict["1_sigma"] = np.percentile(value_array, [15.89, 84.1], axis=0)

    # result_dict["2_sigma"] = np.percentile(value_array, [2.27, 97.725], axis=0)

    return result_dict


def get_histogram_data(bins, data_array, weight_array):
    """
    Function to get the histogram data.

    Also returns the truncated bins where the ends containin only zeros are chopped off
    """

    # Determine the mass bins
    bin_size = np.diff(bins)
    bincenter = (bins[1:] + bins[:-1]) / 2

    # bin and take into account the divison by mass
    hist = np.histogram(data_array, bins=bins, weights=weight_array)[0]

    # Select the non-zero bins and split off the empty ones
    # NOTE: without this, toms method does not work
    non_zero_bins_indices = np.nonzero(hist)[0]

    if non_zero_bins_indices.size != 0:
        truncated_bins = bins[
            non_zero_bins_indices.min() : non_zero_bins_indices.max() + 1
        ]
    else:
        truncated_bins = bins
    return hist, bincenter, truncated_bins


def run_bootstrap(bins, bin_centers, rates, masses, bootstraps=50, verbose=False):
    """
    Function to multiprocess the bootstrapping
    """

    # Get a list of indices
    indices = np.arange(len(rates))

    #########
    # Set up bootstrap array for rates:
    bootstrapped_hist_vals = np.zeros((bootstraps, len(bin_centers)))  # center_bins

    ##########
    # Run bootstrap loop
    for bootstrap_i in range(bootstraps):
        if verbose:
            print("Bootstrap {}".format(bootstrap_i))

        ##############################
        # Get bootstrap indices
        boot_index = np.random.choice(
            indices,
            size=len(indices),
            replace=True,
        )

        #########
        # Calculate rates data with the bootstrapped set of indices
        # Select the quantity values with these indices
        bootstrapped_masses = masses[boot_index]

        # Select the rate values with these indices
        bootstrapped_rates = rates[boot_index]

        ##############################
        # Calculate the rate histogram
        (
            bootstrapped_hist,
            _,
            _,
        ) = get_histogram_data(
            bins=bins,
            data_array=bootstrapped_masses,
            weight_array=bootstrapped_rates,
        )

        # Store unfiltered rate in array
        bootstrapped_hist_vals[bootstrap_i] = bootstrapped_hist

    ###########
    # Calculate median and percentiles
    bootstrapped_median_percentiles_dict = get_median_percentiles(
        bootstrapped_hist_vals
    )

    return bootstrapped_median_percentiles_dict


def plot_bootstrapped_data(
    fig,
    ax,
    bin_centers,
    bin_edges,
    median_percentile_data,
    label,
    color_i="black",
    linestyle_i="solid",
    include_hist_step=False
):
    """ """

    linewidth = 3
    hist_step_alpha = 0.5

    # Plot median and bootstrap
    ax.plot(
        bin_centers,
        median_percentile_data["median"][0],
        lw=linewidth,
        c=color_i,
        zorder=13,
        linestyle=linestyle_i,
        label=label,
    )

    # fill between for the bounds
    ax.fill_between(
        bin_centers,
        median_percentile_data["90%_CI"][0],
        median_percentile_data["90%_CI"][1],
        alpha=0.4,
        zorder=11,
        color=color_i,
    )  # 1-sigma

    # Plot step histogram
    if include_hist_step:
        ax.hist(
            bin_centers,
            weights=median_percentile_data["median"][0],
            bins=bin_edges,
            histtype="step",
            lw=linewidth,
            color=color_i,
            zorder=200,
            alpha=hist_step_alpha,
            linestyle=linestyle_i,
        )
