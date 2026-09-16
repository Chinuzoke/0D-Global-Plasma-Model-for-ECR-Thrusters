"""
species.py
----------
Propellant / plasma species properties.

Both reference models (Yamamoto-group gridded ECR ion thrusters: mu1, mu10,
MIPS/exospheric thrusters; Jorns-group PEPL magnetic-nozzle ECR thruster,
Wachs & Jorns, AIAA-2021-3382) are run almost exclusively on xenon, so Xe is
the default species. Krypton and Argon are included for completeness since
some Yamamoto-group papers (e.g. water/H2/He microplasma thruster work) and
general ECR/global-model literature explore alternative propellants.

Ionization and (single, lumped) excitation energies are taken from standard
NIST atomic spectra data and from the values used inside Lieberman &
Lichtenberg-style global models (the common numerical basis cited directly
in Wachs & Jorns 2021, Ref. [10]-[11] of that paper).
"""

from dataclasses import dataclass
from constants import AMU


@dataclass(frozen=True)
class Species:
    name: str
    mass_amu: float          # atomic mass [amu]
    E_iz_eV: float           # ionization energy (ground state -> singly ionized) [eV]
    E_ex_eV: float           # representative (lumped) excitation energy [eV]

    @property
    def mass_kg(self):
        return self.mass_amu * AMU


# ---- Propellant database -----------------------------------------------
XENON = Species(name="Xe", mass_amu=131.293, E_iz_eV=12.13, E_ex_eV=11.6)
KRYPTON = Species(name="Kr", mass_amu=83.798, E_iz_eV=14.00, E_ex_eV=10.6)
ARGON = Species(name="Ar", mass_amu=39.948, E_iz_eV=15.76, E_ex_eV=11.55)

# Electron "species" is needed for mass ratio calculations (e.g. sheath
# potential Vs = (Te/2) * ln(mi / (2*pi*me)) used identically in both
# reference models).
from constants import M_ELECTRON  # noqa: E402

ELECTRON_MASS_KG = M_ELECTRON

PROPELLANTS = {
    "Xe": XENON,
    "Kr": KRYPTON,
    "Ar": ARGON,
}


def get_species(name: str) -> Species:
    try:
        return PROPELLANTS[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown propellant '{name}'. Available: {list(PROPELLANTS)}"
        ) from exc
