import copy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from syntheticstellarpopconvolve.general_functions import (  # calculate_bin_edges,
    calculate_bincenters,
)
from syntheticstellarpopconvolve.metallicity_distributions import (
    metallicity_distribution_vanSon2022,
)
from syntheticstellarpopconvolve.SFR_dict_plotting_routines import plot_sfr_dict
from syntheticstellarpopconvolve.starformation_rate_distributions import (
    starformation_rate_distribution_vanSon2023,
)

# redshift_bin_edges = np.arange(0, 8, 0.05)
# log_metallicity_bin_edges = np.arange(-8, 0, 0.05)

# Set up redshift bin info
num_redshift_bins = 100
redshift_bin_edges = np.linspace(0, 8, num_redshift_bins)
redshift_bin_centers = calculate_bincenters(redshift_bin_edges)

# Set up metallicity bin info
num_metallicity_bins = 200
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
    "metallicity_distribution_array": dpdlogZ,
}

axis_dict = plot_sfr_dict(
    high_res_sfr_dict,
    time_type="redshift",
    metallicity_string="logZ",
    metallicity_distribution_multiply_by_metallicity_bin_sizes=True,
    metallicity_distribution_multiply_by_sfr=False,
    metallicity_distribution_scale="log10",
    metallicity_distribution_cmap=copy.copy(plt.cm.viridis),
    return_axis_dict=True,
    figsize=(8, 8),
    fontsize=12,
)

plt.show()
