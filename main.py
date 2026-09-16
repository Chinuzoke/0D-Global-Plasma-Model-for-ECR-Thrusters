
"""
main.py
-------
Main entry point. Runs BOTH cross-referenced 0D global-model architectures
to steady state, prints plasma properties and performance metrics, and
generates plots. 

Case A: Jorns/Wachs magnetic-nozzle ECR thruster
    Operating point: 30 W absorbed power, 1 sccm xenon -- the baseline
    condition used throughout Wachs & Jorns, AIAA-2021-3382.
    Geometry (R, L) is NOT given numerically in the excerpted text of that
    paper (it is a small coaxial ONERA-style ECR thruster, similar to
    Cannat et al. 2015 / Vialis et al. 2018); we use illustrative
    small-thruster dimensions (R=1.5 cm, L=3 cm) typical of that thruster
    class and flag this clearly -- swap in exact dimensions if you have
    them for a publication-grade reproduction.

Case B: Yamamoto-lineage gridded ECR ion thruster (mu10-class)
    Operating point: illustrative low-power mu10-class condition
    (~30-35 W discharge, small sccm-level xenon flow, ~1500 V beam
    voltage), consistent with the general performance envelope reported
    across the Yamamoto/Tsukizaki/Nishiyama mu-series papers referenced in
    the module docstrings. Exact ISAS proprietary chamber dimensions are
    not published in the fetched sources; illustrative dimensions
    (R=5 cm, L=7 cm, screen open-area fraction 0.67, typical of a
    multi-aperture screen grid) are used and clearly flagged.

Both cases are examples: change the parameter blocks below (or import
this package's modules directly) to model your own thruster geometry and
operating conditions.
"""

import numpy as np

from geometry import ThrusterGeometry
from species import XENON
from cross_sections import RateCoefficientModel
from models.jorns_magnetic_nozzle import JornsParams, rhs as jorns_rhs, performance as jorns_performance
from models.yamamoto_gridded import YamamotoParams, rhs as yam_rhs, performance as yam_performance
from solver import run_to_steady_state
import plotting


def run_case_jorns():
    print("=" * 70)
    print("CASE A: Jorns/Wachs magnetic-nozzle ECR thruster (AIAA-2021-3382)")
    print("=" * 70)

    geom = ThrusterGeometry(R=0.015, L=0.03, screen_open_area_frac=1.0)
    rate_model = RateCoefficientModel(XENON)

    params = JornsParams(
        geometry=geom,
        species=XENON,
        rate_model=rate_model,
        mdot_sccm=1.0,
        Tg_K=500.0,
        power_profile=lambda t: 30.0,   # constant 30 W CW, baseline case
        h_R=0.5,
        h_B=0.4,
        h_T=0.5,
        M_det=2.0,
        V_eff_frac=1.0,
    )

    # Sensible initial guess: modest ionization fraction, cold-ish gas fill
    ng0 = params.mdot_particles_per_s / (params.geometry.area_end
                                         * np.sqrt(8 * 1.380649e-23 * params.Tg_K
                                                   / (np.pi * params.mi)) / 4.0)
    y0 = [1e16, ng0, 3.0]

    sol, y_ss, converged = run_to_steady_state(
        jorns_rhs, y0, params, t_span=(0.0, 2e-3))
    perf = jorns_performance(y_ss, params)

    print("  Performance:")
    for k, v in perf.items():
        print(f"    {k:20s} = {v:.4g}")

    return sol, y_ss, params, perf


def run_case_yamamoto():
    print("=" * 70)
    print("CASE B: Yamamoto-lineage gridded ECR ion thruster (mu10-class)")
    print("=" * 70)

    geom = ThrusterGeometry(R=0.05, L=0.07, screen_open_area_frac=0.67)
    rate_model = RateCoefficientModel(XENON)

    params = YamamotoParams(
        geometry=geom,
        species=XENON,
        rate_model=rate_model,
        mdot_sccm=1.5,
        Tg_K=500.0,
        power_profile=lambda t: 35.0,
        h_R=0.3,
        h_B=0.4,
        h_S=0.6,
        V_eff_frac=1.0,
        V_beam=1500.0,
    )

    ng0 = params.mdot_particles_per_s / (params.geometry.area_screen_open
                                         * np.sqrt(8 * 1.380649e-23 * params.Tg_K
                                                   / (np.pi * params.mi)) / 4.0)
    y0 = [1e16, ng0, 5.0]

    sol, y_ss, converged = run_to_steady_state(
        yam_rhs, y0, params, t_span=(0.0, 2e-3))
    perf = yam_performance(y_ss, params)

    print("  Performance:")
    for k, v in perf.items():
        print(f"    {k:20s} = {v:.4g}")

    return sol, y_ss, params, perf


if __name__ == "__main__":
    sol_j, yss_j, params_j, perf_j = run_case_jorns()
    print()
    sol_y, yss_y, params_y, perf_y = run_case_yamamoto()

    print()
    print("Generating plots in ./outputs_plots/ ...")
    plotting.plot_transient(sol_j, title="Jorns/Wachs ECR thruster - transient approach to steady state",
                            filename="jorns_transient.png")
    plotting.plot_transient(sol_y, title="Yamamoto-lineage gridded ECR thruster - transient approach to steady state",
                            filename="yamamoto_transient.png")
    plotting.plot_hR_sweep(params_j, filename="jorns_hR_sweep.png")
    print("Done.")
