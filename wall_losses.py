"""
wall_losses.py
--------------
Sheath-edge physics and edge-to-center density ("h factor") machinery shared
by both reference 0D global models.

Directly implements the quantities appearing in Wachs & Jorns
(AIAA-2021-3382), Sec. II, Eqs. (1)-(7):

    u_B          = sqrt(e Te / m_i)                     Bohm speed
    V_s          = (Te/2) * ln( m_i / (2 pi m_e) )       sheath potential [V]
    v_therm      = sqrt(8 kB Tg / (pi m_i))              neutral thermal speed
    h_R, h_L, h_T                                        edge-to-center
                                                          density ratios for
                                                          the radial wall,
                                                          back wall, and
                                                          exit/throat plane

as well as the classical Godyak/Lieberman-Lichtenberg low-pressure-limit
closed forms for h_l and h_R as a function of the ion-neutral mean free
path (Lieberman & Lichtenberg 2005, Eqs. 5.53-5.55, 10.3-10.7), which is
the formula set the Yamamoto-group gridded-thruster lineage and most
general ECR/helicon/ICP global models (e.g. the Lafleur helicon model cited
by Wachs & Jorns as Ref. [7]) use to CLOSE h_R/h_L when they are not simply
left as free/fitted parameters.

Wachs & Jorns explicitly treat h_R as a free/swept parameter standing in for
magnetic-field confinement physics not captured by the neutral-collisional
Godyak picture (their Sec. III.A, Fig. 2a). This module supports BOTH usage
modes:
  - `godyak_h_factors(...)`  -> analytic collisional closure (Yamamoto-style
     / general low-pressure global-model closure)
  - directly passing a constant h_R (Jorns/Wachs-style free parameter, e.g.
     representative magnetically-confined values h_R ~ 0.2-0.5)
"""

import numpy as np
from constants import E_CHARGE, M_ELECTRON, K_BOLTZMANN, eV_to_J


def bohm_velocity(Te_eV, m_i_kg):
    """u_B = sqrt(e*Te / m_i)  [m/s]   (Te in eV)."""
    return np.sqrt(eV_to_J(Te_eV) / m_i_kg)


def sheath_potential_V(Te_eV, m_i_kg):
    """
    Floating-sheath potential drop V_s [in units of Te, i.e. already
    multiplied through] following Wachs & Jorns Eq. (6)/(7):

        V_s = (Te / 2) * ln( m_i / (2 pi m_e) )      [V, if Te in eV]

    This is the standard planar-sheath floating potential result (see also
    Lieberman & Lichtenberg 2005, Eq. 2.4.14).
    """
    return 0.5 * Te_eV * np.log(m_i_kg / (2.0 * np.pi * M_ELECTRON))


def neutral_thermal_velocity(Tg_K, m_i_kg):
    """v_therm = sqrt(8 kB Tg / (pi m_i))  [m/s], Wachs & Jorns Eq. (2)."""
    return np.sqrt(8.0 * K_BOLTZMANN * Tg_K / (np.pi * m_i_kg))


def godyak_h_factors(R, L, lambda_i):
    """
    Low-pressure-limit collisional edge-to-center density ratios
    (Lieberman & Lichtenberg 2005, Sec. 5.3 / 10.2; the standard closure
    used across the ECR/helicon/ICP thruster global-model literature,
    e.g. by the Yamamoto-group lineage and by Lafleur's helicon model that
    Wachs & Jorns cite as their Ref. [7]).

    Parameters
    ----------
    R, L : float
        Chamber radius and length [m].
    lambda_i : float
        Ion-neutral (charge-exchange/momentum-transfer) mean free path [m].

    Returns
    -------
    h_L, h_R : float
        Axial (end-wall) and radial (side-wall) edge-to-center density
        ratios (dimensionless, 0 < h < 1).
    """
    h_L = 0.86 / np.sqrt(3.0 + (L / (2.0 * lambda_i)))
    h_R = 0.80 / np.sqrt(4.0 + (R / lambda_i))
    return h_L, h_R


def ion_neutral_mfp(n_g, sigma_cx=1.0e-18):
    """
    Ion-neutral charge-exchange mean free path lambda_i = 1/(n_g * sigma_cx).

    sigma_cx defaults to a representative xenon ion-neutral charge-exchange
    cross section (~1e-18 m^2 order of magnitude at thermal-ish relative
    energies; consistent with values compiled across the Hall/ion-thruster
    global-model literature). Override with a tabulated/energy-dependent
    value for higher fidelity.
    """
    return 1.0 / (n_g * sigma_cx)
