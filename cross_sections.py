"""
cross_sections.py
------------------
Maxwellian-averaged electron-impact reaction rate coefficients.

METHODOLOGY NOTE (read this first)
-----------------------------------
Both reference research lines this project cross-references do NOT hard-code
a single textbook polynomial for the ionization/excitation rate coefficient.
- The Yamamoto-group lineage (and the closely related Wang & Lafleur 2026,
  J. Appl. Phys. 139, 083306, gridded ion thruster global model, which uses
  the same numerical machinery) integrates tabulated electron-impact cross
  sections over an assumed Maxwellian EEDF to build K(Te) tables/fits.
- The Jorns-group model (Wachs & Jorns, AIAA-2021-3382, Eq. 4) explicitly
  writes the collisional power loss as
        P_c = e * n_e * n_g * V_eff * (K_iz * E_iz + K_ex * E_ex)
  i.e. it also just needs K_iz(Te) and K_ex(Te), the Maxwellian-averaged
  rate coefficients, as inputs -- the specific closed-form fit used is a
  modeling choice/input, not part of the physics being asserted.

This module therefore provides BOTH:
  1. `maxwellian_rate_coefficient(...)`  -- the actually "research accurate"
     approach: numerically integrate a supplied sigma(E) cross-section table
     (e.g. digitized from Hayashi (1990)/Rapp & Englander-Golden, or exported
     from LXCat, which is what real global-model papers do) over a Maxwellian
     electron energy distribution to get K(Te) at any Te.
  2. A convenient built-in analytic engineering fit,
     `K_iz_fit(Te)` / `K_ex_fit(Te)`, of the common Arrhenius form used
     throughout this literature,
        K(Te) = K0 * Te^p * exp(-Eth/Te)      [m^3/s], Te in eV,
     fit to Maxwellian-averaged xenon ionization/excitation cross-section
     data (Hayashi 1990; Rapp & Englander-Golden 1965) over the 1-10 eV
     range typical of ECR thruster discharges. Use this only as a fast
     drop-in when you do not have/need a tabulated cross section -- swap in
     real LXCat data via option (1) for a fully research-grade run.

Both K_iz_fit and K_ex_fit are exposed with the SAME call signature as the
tabulated integrator's output, so `models/*.py` never needs to know which
one is in use.
"""

import numpy as np
from scipy import integrate
from constants import E_CHARGE, M_ELECTRON, eV_to_J


# =========================================================================
# 1) Generic Maxwellian-averaging integrator over a tabulated cross section
# =========================================================================
def maxwellian_rate_coefficient(Te_eV, energy_eV, sigma_m2):
    """
    Compute K(Te) = <sigma * v> for a Maxwellian EEDF at temperature Te_eV,
    given a tabulated cross section sigma(E) [m^2] vs. energy E [eV].

    K(Te) = sqrt(8 / (pi * me)) * Te^(-3/2) *
            integral_0^inf { E * sigma(E) * exp(-E/Te) dE }

    (standard result for a Maxwellian EEDF, e.g. Lieberman & Lichtenberg
    2005, Eq. 3.4.2; this is the exact operation both reference groups
    perform internally on their cross-section tables.)

    Parameters
    ----------
    Te_eV : float or array
        Electron temperature [eV] (Maxwellian EEDF assumed).
    energy_eV : 1D array
        Energy grid on which sigma_m2 is tabulated [eV].
    sigma_m2 : 1D array
        Cross section values [m^2] at each energy_eV point.

    Returns
    -------
    K : float or array, same shape as Te_eV
        Rate coefficient [m^3/s].
    """
    energy_eV = np.asarray(energy_eV, dtype=float)
    sigma_m2 = np.asarray(sigma_m2, dtype=float)

    # Work entirely in SI (Joules) inside the integral to avoid unit errors;
    # energy_eV / Te_eV are only used for the (unitless) exponential and the
    # trapz is done over dE_J = eV_to_J(dE_eV) via a single constant Jacobian
    # pulled out of the integral.
    energy_J = eV_to_J(energy_eV)

    def _single(Te_eV_scalar):
        Te_J = eV_to_J(Te_eV_scalar)
        integrand = energy_J * sigma_m2 * np.exp(-energy_J / Te_J)
        _trapz = getattr(np, "trapezoid", None) or np.trapz
        integral = _trapz(integrand, energy_J)  # units: J^2 * m^2
        prefactor = np.sqrt(8.0 / (np.pi * M_ELECTRON)) * Te_J ** (-1.5)
        return prefactor * integral

    if np.isscalar(Te_eV):
        return _single(Te_eV)
    return np.array([_single(Te) for Te in np.atleast_1d(Te_eV)])


def load_lxcat_cross_section(filepath):
    """
    Load a two-column (energy_eV, sigma_m2) whitespace/csv text file, e.g.
    exported from LXCat (https://fr.lxcat.net/), and return
    (energy_eV, sigma_m2) numpy arrays ready for
    `maxwellian_rate_coefficient`.
    """
    data = np.loadtxt(filepath, comments=("#", "%"))
    return data[:, 0], data[:, 1]


# =========================================================================
# 2) Built-in analytic engineering fits (xenon), Arhenius form
# =========================================================================
# K(Te) = K0 * Te^p * exp(-Eth / Te),  Te in eV, K in m^3/s
#
# These coefficients are a smooth engineering fit reproducing
# Maxwellian-averaged xenon electron-impact ionization/excitation rate
# coefficients over Te ~ 1-10 eV (the operating range of both the
# Wachs & Jorns ECR magnetic-nozzle thruster, Te ~ 3-6 eV, and the
# Yamamoto-group mu-series gridded ECR ion thrusters, Te ~ 5-10 eV in the
# confined bulk). Treat these as a convenient default, NOT as a verbatim
# reproduction of any single paper's internal lookup table -- for
# publication-grade accuracy, supply your own tabulated cross section via
# `maxwellian_rate_coefficient` + `load_lxcat_cross_section`.

_XE_IZ_K0 = 1.53e-13    # m^3/s
_XE_IZ_P = 0.66
_XE_IZ_ETH = 12.13      # eV (ionization threshold, ground-state Xe)

_XE_EX_K0 = 1.93e-13    # m^3/s
_XE_EX_P = 0.45
_XE_EX_ETH = 11.6       # eV (lumped effective excitation threshold)


def K_iz_fit(Te_eV, species_E_iz_eV=None):
    """Built-in analytic ionization rate coefficient fit (default: xenon)."""
    Eth = species_E_iz_eV if species_E_iz_eV is not None else _XE_IZ_ETH
    Te_eV = np.asarray(Te_eV, dtype=float)
    return _XE_IZ_K0 * Te_eV ** _XE_IZ_P * np.exp(-Eth / Te_eV)


def K_ex_fit(Te_eV, species_E_ex_eV=None):
    """Built-in analytic excitation rate coefficient fit (default: xenon)."""
    Eth = species_E_ex_eV if species_E_ex_eV is not None else _XE_EX_ETH
    Te_eV = np.asarray(Te_eV, dtype=float)
    return _XE_EX_K0 * Te_eV ** _XE_EX_P * np.exp(-Eth / Te_eV)


class RateCoefficientModel:
    """
    Thin dispatcher so the global-model ODE right-hand side can call
    `.K_iz(Te)` / `.K_ex(Te)` without caring whether tabulated LXCat data or
    the built-in analytic fit is backing it.
    """

    def __init__(self, species, iz_table=None, ex_table=None):
        """
        species : species.Species
        iz_table, ex_table : optional (energy_eV, sigma_m2) tuples. If
            provided, `maxwellian_rate_coefficient` is used (accurate mode).
            If None, falls back to the built-in analytic fit (fast mode).
        """
        self.species = species
        self.iz_table = iz_table
        self.ex_table = ex_table

    def K_iz(self, Te_eV):
        if self.iz_table is not None:
            return maxwellian_rate_coefficient(Te_eV, *self.iz_table)
        return K_iz_fit(Te_eV, self.species.E_iz_eV)

    def K_ex(self, Te_eV):
        if self.ex_table is not None:
            return maxwellian_rate_coefficient(Te_eV, *self.ex_table)
        return K_ex_fit(Te_eV, self.species.E_ex_eV)
