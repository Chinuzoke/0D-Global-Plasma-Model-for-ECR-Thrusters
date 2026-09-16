"""
constants.py
------------
Fundamental physical constants used throughout the 0D global model.

All quantities are SI unless explicitly noted (electron temperature Te is
carried in eV throughout this codebase, as is standard in the plasma
propulsion literature -- see Lieberman & Lichtenberg, "Principles of Plasma
Discharges and Materials Processing", 2nd ed., Wiley, 2005, Ch. 10, which is
the common numerical ancestor of both the Yamamoto-group gridded ECR ion
thruster models and the Jorns-group (PEPL) magnetic-nozzle ECR thruster
model implemented here.
"""

# ---- Universal constants (SI) ----------------------------------------
E_CHARGE = 1.602176634e-19      # Elementary charge [C]
M_ELECTRON = 9.1093837015e-31   # Electron mass [kg]
K_BOLTZMANN = 1.380649e-23      # Boltzmann constant [J/K]
EPS_0 = 8.8541878128e-12        # Vacuum permittivity [F/m]
AMU = 1.66053906660e-27         # Atomic mass unit [kg]
C_LIGHT = 2.99792458e8          # Speed of light [m/s]

# ---- Convenience conversions -------------------------------------------
def eV_to_J(x_eV):
    """Convert an energy in eV to Joules."""
    return x_eV * E_CHARGE


def J_to_eV(x_J):
    """Convert an energy in Joules to eV."""
    return x_J / E_CHARGE


def Te_eV_to_K(Te_eV):
    """Convert an electron temperature in eV to Kelvin (T = E/k_B)."""
    return Te_eV * E_CHARGE / K_BOLTZMANN


# ---- Standard gas conditions -------------------------------------------
SCCM_TO_PARTICLES_PER_S = 4.478e17
# 1 standard cm^3/min of gas = 4.478e17 particles/s (species-independent;
# this is the standard SCCM definition used in both the Wachs & Jorns 2021
# AIAA-2021-3382 paper and the Yamamoto-group ion-thruster papers when
# converting mass flow rate mdot = (sccm * SCCM_TO_PARTICLES_PER_S) * m_i).

TORR_TO_PA = 133.322
