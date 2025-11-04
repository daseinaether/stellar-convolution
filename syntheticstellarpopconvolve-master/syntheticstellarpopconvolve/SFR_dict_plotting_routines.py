import copy
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

from syntheticstellarpopconvolve.check_and_update_sfr_dict import (
    check_and_update_sfr_dict,
)


def load_mpl_rc():  # DH0001

    # https://matplotlib.org/users/customizing.html
    mpl.rc(
        "axes",
        labelweight="normal",
        linewidth=2,
        labelsize=30,
        grid=True,
        titlesize=40,
        facecolor="white",
    )

    mpl.rc("savefig", dpi=100)

    mpl.rc("lines", linewidth=4, color="g", markeredgewidth=2)

    mpl.rc(
        "ytick",
        **{
            "labelsize": 30,
            "color": "k",
            "left": True,
            "right": True,
            "major.size": 12,
            "major.width": 2,
            "minor.size": 6,
            "minor.width": 2,
            "major.pad": 12,
            "minor.visible": True,
            "direction": "inout",
        },
    )

    mpl.rc(
        "xtick",
        **{
            "labelsize": 30,
            "top": True,
            "bottom": True,
            "major.size": 12,
            "major.width": 2,
            "minor.size": 6,
            "minor.width": 2,
            "major.pad": 12,
            "minor.visible": True,
            "direction": "inout",
        },
    )

    mpl.rc("legend", frameon=False, fontsize=30, title_fontsize=30)

    mpl.rc("contour", negative_linestyle="solid")

    mpl.rc(
        "figure",
        figsize=[16, 16],
        titlesize=30,
        dpi=100,
        facecolor="white",
        edgecolor="white",
        frameon=True,
        max_open_warning=10,
        # autolayout=True
    )

    mpl.rc(
        "legend",
        fontsize=20,
        handlelength=2,
        loc="best",
        fancybox=False,
        numpoints=2,
        framealpha=None,
        scatterpoints=3,
        edgecolor="inherit",
    )

    mpl.rc("savefig", dpi="figure", facecolor="white", edgecolor="white")

    mpl.rc("grid", color="b0b0b0", alpha=0.5)

    mpl.rc("image", cmap="viridis")

    mpl.rc(
        "font",
        # weight='bold',
        serif="Palatino",
        size=20,
    )

    mpl.rc("errorbar", capsize=2)

    mpl.rc("mathtext", default="sf")


def plot_sfr_dict(  # DH0001
    sfr_dict,
    time_type,
    metallicity_string="Z",
    metallicity_distribution_scale="linear",
    metallicity_distribution_multiply_by_metallicity_bin_sizes=False,
    metallicity_distribution_multiply_by_sfr=False,
    metallicity_distribution_max_logdiff=3,
    metallicity_distribution_cmap=None,
    return_axis_dict=False,
    figsize=None,
    fontsize=20,
):
    """
    Function to plot the star formation rate
    """

    load_mpl_rc()

    # check if there is any metallicity info
    has_metallicity_info = (
        "metallicity_bin_edges" in sfr_dict
        and "metallicity_distribution_array" in sfr_dict
    )

    # check and update the sfr dict
    from astropy.cosmology import Planck13 as cosmo  # Planck 2013

    sfr_dict = check_and_update_sfr_dict(
        sfr_dict=sfr_dict,
        config={
            "logger": logging.getLogger(__name__),
            "time_type": time_type,
            "cosmology": cosmo,
        },
        requires_name=False,
        requires_metallicity_info=has_metallicity_info,
        time_type=time_type,
    )

    ########
    # Set up figure canvas
    axis_dict = {}
    if figsize is not None:
        axis_dict["fig"] = plt.figure(figsize=figsize)
    else:
        axis_dict["fig"] = plt.figure(figsize=(20, 16 if has_metallicity_info else 6))

    if has_metallicity_info:
        gs = axis_dict["fig"].add_gridspec(nrows=3, ncols=8)

        axis_dict["ax_sfr"] = axis_dict["fig"].add_subplot(gs[0, :-2])
        axis_dict["ax_mssfr"] = axis_dict["fig"].add_subplot(gs[1:, :-2])
        axis_dict["ax_bar"] = axis_dict["fig"].add_subplot(gs[1:, -1])

        axis_dict["ax_sfr"].tick_params(axis="both", which="major", labelsize=fontsize)
        axis_dict["ax_mssfr"].tick_params(
            axis="both", which="major", labelsize=fontsize
        )
        axis_dict["ax_bar"].tick_params(axis="both", which="major", labelsize=fontsize)

    else:
        gs = axis_dict["fig"].add_gridspec(nrows=1, ncols=2)

        axis_dict["ax_sfr"] = axis_dict["fig"].add_subplot(gs[0, :])
        axis_dict["ax_sfr"].tick_params(axis="both", which="major", labelsize=fontsize)

    ########
    # Plot Star Formation Rate (SFR)

    if "time_bin_edges" not in sfr_dict:
        raise ValueError("'time_bin_edges' are not present in the sfr dict!")
    if "starformation_rate_array" not in sfr_dict:
        raise ValueError("'starformation_rate_array' is not present in the sfr dict!")

    time_bin_edges = sfr_dict["time_bin_centers"]
    sfr_array = sfr_dict["starformation_rate_array"]

    axis_dict["ax_sfr"].plot(time_bin_edges, sfr_array)

    if time_type == "lookback_time":
        axis_dict["ax_sfr"].set_xlabel(
            "Lookback Time [{}]".format(
                sfr_dict["time_bin_centers"].unit.to_string("latex_inline")
            ),
            fontsize=fontsize,
        )
    if time_type == "redshift":
        axis_dict["ax_sfr"].set_xlabel("Redshift", fontsize=fontsize)

    axis_dict["ax_sfr"].set_ylabel(
        "Star Formation Rate\n[{}]".format(
            sfr_dict["starformation_rate_array"].unit.to_string("latex_inline")
        ),
        fontsize=fontsize,
    )
    axis_dict["ax_sfr"].set_title("Star Formation Rate", fontsize=fontsize)

    # set order of magnitude text lower
    text = axis_dict["ax_sfr"].yaxis.get_offset_text()  # Get the text object
    text.set_size(fontsize)  # # Set the size.

    ########
    # Plot the metallicity
    if has_metallicity_info:
        # Set up the meshgrid that we will use
        time_mesh, metallicity_mesh = np.meshgrid(
            # sfr_dict["time_bin_centers"], sfr_dict["metallicity_bin_centers"]
            sfr_dict["time_bin_edges"],
            sfr_dict["metallicity_bin_edges"],
        )

        cmap_metallicity_fraction_hist = (
            copy.copy(plt.cm.jet)
            if metallicity_distribution_cmap is None
            else metallicity_distribution_cmap
        )

        #
        z_vals = sfr_dict["metallicity_distribution_array"]
        if metallicity_distribution_multiply_by_metallicity_bin_sizes:
            z_vals = z_vals * sfr_dict["metallicity_bin_sizes"][np.newaxis, :]
        if metallicity_distribution_multiply_by_sfr:
            z_vals = z_vals * sfr_array[:, np.newaxis].value
        vmin = z_vals[~np.isnan(z_vals)].min()
        vmax = z_vals[~np.isnan(z_vals)].max()

        if metallicity_distribution_scale == "linear":
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
        elif metallicity_distribution_scale == "log10":
            # DIFF_MAX_THRESHOLD = 1
            DIFF_MAX_THRESHOLD = metallicity_distribution_max_logdiff

            vmin = 10 ** (np.log10(vmax) - DIFF_MAX_THRESHOLD)
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)

        if time_type == "lookback_time":
            _ = axis_dict["ax_mssfr"].pcolormesh(
                time_mesh.value,
                metallicity_mesh,
                z_vals.T,
                norm=norm,
                cmap=cmap_metallicity_fraction_hist,
                # shading="auto",
                antialiased=True,
                rasterized=True,
                shading="flat",
            )

        if time_type == "redshift":
            _ = axis_dict["ax_mssfr"].pcolormesh(
                time_mesh,
                metallicity_mesh,
                z_vals.T,
                norm=norm,
                cmap=cmap_metallicity_fraction_hist,
                # shading="auto",
                antialiased=True,
                rasterized=True,
                shading="flat",
            )

        # make colorbar
        _ = mpl.colorbar.ColorbarBase(
            axis_dict["ax_bar"], norm=norm, cmap=cmap_metallicity_fraction_hist
        )
        # cbar.ax.set_ylabel(r"Fraction of ZAMS mass")
        # axis_dict['fig'].colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_metallicity_fraction_hist), ax=axis_dict['ax_bar'])

        #
        axis_dict["ax_sfr"].set_xlabel("")
        axis_dict["ax_sfr"].set_xticklabels([])
        if time_type == "lookback_time":
            axis_dict["ax_mssfr"].set_xlabel(
                "Lookback Time [{}]".format(
                    sfr_dict["time_bin_centers"].unit.to_string("latex_inline")
                ),
                fontsize=fontsize,
            )
        if time_type == "redshift":
            axis_dict["ax_mssfr"].set_xlabel("Redshift", fontsize=fontsize)

        ax_mssfr_title = r"$dP/d{}$".format(metallicity_string)

        if (
            metallicity_distribution_multiply_by_metallicity_bin_sizes
            or metallicity_distribution_multiply_by_sfr
        ):
            ax_mssfr_title = "({})".format(ax_mssfr_title)
        if metallicity_distribution_multiply_by_metallicity_bin_sizes:
            ax_mssfr_title = "{}{}".format(
                ax_mssfr_title, "*Δ{}".format(metallicity_string)
            )
        if metallicity_distribution_multiply_by_sfr:
            ax_mssfr_title = "{}{}".format(ax_mssfr_title, "*SFR")

        axis_dict["ax_mssfr"].set_title(ax_mssfr_title, fontsize=fontsize)

        axis_dict["ax_mssfr"].set_ylabel(
            "Metallicity [{}]".format(metallicity_string), fontsize=fontsize
        )
        axis_dict["ax_mssfr"].set_ylim(
            [
                sfr_dict["metallicity_bin_edges"].min(),
                sfr_dict["metallicity_bin_edges"].max(),
            ]
        )

        #
        axis_dict["fig"].align_ylabels([axis_dict["ax_sfr"], axis_dict["ax_mssfr"]])

    #############
    #
    if return_axis_dict:
        return axis_dict
    else:
        plt.show()


if __name__ == "__main__":

    import astropy.units as u

    from syntheticstellarpopconvolve.general_functions import calculate_bincenters
    from syntheticstellarpopconvolve.metallicity_distributions import (
        metallicity_distribution_vanSon2022,
    )
    from syntheticstellarpopconvolve.starformation_rate_distributions import (
        starformation_rate_distribution_vanSon2023,
    )

    # Set up redshift bin info
    num_redshift_bins = 400
    redshift_bin_edges = np.linspace(0, 8, num_redshift_bins)
    redshift_bin_centers = calculate_bincenters(redshift_bin_edges)

    # Set up metallicity bin info
    num_metallicity_bins = 500
    log_metallicity_bin_edges = np.linspace(-12, 0, num_metallicity_bins)
    log_metallicity_bin_centers = calculate_bincenters(log_metallicity_bin_edges)

    #
    sfr = starformation_rate_distribution_vanSon2023(redshift_bin_centers).to(
        u.Msun / u.yr / u.Gpc**3
    )

    #
    dpdlogZ = metallicity_distribution_vanSon2022(
        log_metallicity_centers=log_metallicity_bin_centers,
        redshifts=redshift_bin_centers,
    )

    high_res_sfr_dict = {
        "redshift_bin_edges": redshift_bin_edges,
        "starformation_rate_array": sfr,
        "metallicity_bin_edges": log_metallicity_bin_edges,
        "metallicity_distribution_array": dpdlogZ,  # We need to transpose!
    }

    axis_dict = plot_sfr_dict(
        high_res_sfr_dict,
        time_type="redshift",
        metallicity_string="logZ",
        metallicity_distribution_multiply_by_metallicity_bin_sizes=True,
        metallicity_distribution_multiply_by_sfr=True,
        metallicity_distribution_scale="log10",
        metallicity_distribution_cmap=copy.copy(plt.cm.viridis),
        return_axis_dict=True,
        figsize=(8, 8),
        fontsize=12,
    )
    plt.show()
