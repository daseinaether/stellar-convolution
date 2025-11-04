"""
Functions to calculate redshift shell info
"""

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve.cosmology_utils import redshift_to_lookback_time
from syntheticstellarpopconvolve.general_functions import calculate_bincenters


def create_shell_volume_dict(redshift_bin_edges, config):
    """
    Function that can generate a dictionary of shell volumes based on an input center redshift array
    """

    #
    config["logger"].debug("Storing redshift-shell info in SFR dict")

    #
    redshift_bin_centers = calculate_bincenters(redshift_bin_edges)

    # Calculate comoving volumes
    comoving_volumes_at_redshift_edges = (
        config["cosmology"].comoving_volume(redshift_bin_edges).to(u.Gpc**3)
    )

    # Calculate the comoving shell volumes
    comoving_shell_volumes = np.diff(comoving_volumes_at_redshift_edges)

    # create dict with the value at center of 'shell'
    comoving_shell_volumes_dict = {
        redshift_bin_centers[i]: {
            "shell_volume": comoving_shell_volumes[i],
            "lower_edge_shell_redshift": redshift_bin_edges[i],
            "upper_edge_shell_redshift": redshift_bin_edges[i + 1],
            "lower_edge_shell_lookback_time": redshift_to_lookback_time(
                redshift=redshift_bin_edges[i], cosmology=config["cosmology"]
            ),
            "upper_edge_shell_lookback_time": redshift_to_lookback_time(
                redshift=redshift_bin_edges[i + 1], cosmology=config["cosmology"]
            ),
            "delta_shell_redshift": np.abs(
                redshift_bin_edges[i + 1] - redshift_bin_edges[i]
            ),
            "delta_shell_lookback_time": np.abs(
                redshift_to_lookback_time(
                    redshift=redshift_bin_edges[i + 1], cosmology=config["cosmology"]
                )
                - redshift_to_lookback_time(
                    redshift=redshift_bin_edges[i], cosmology=config["cosmology"]
                )
            ),
            "center_shell": redshift_bin_centers[i],
        }
        for i in range(len(redshift_bin_centers))
    }

    return comoving_shell_volumes_dict


def store_redshift_shell_info(config, sfr_dict):
    """
    Function to add the redshift shell info dict to the hdf5 file
    """

    if config["time_type"] == "redshift":
        ##################
        # Create shell volume dict
        redshift_shell_volume_dict = create_shell_volume_dict(
            redshift_bin_edges=sfr_dict["redshift_bin_edges"],
            config=config,
        )
        sfr_dict["redshift_shell_volume_dict"] = redshift_shell_volume_dict

    return sfr_dict
