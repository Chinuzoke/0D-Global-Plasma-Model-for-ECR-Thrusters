"""
plotting.py
-----------
Plotting helpers. Produces:
  - transient time-history plots of n_e, n_g, Te (mirrors Wachs & Jorns
    Fig. 1: "Time resolved global model showing (a) absorbed power
    (b) electron temperature & electron density.")
  - a parameter sweep of efficiency vs. h_R (mirrors Wachs & Jorns Fig. 2a:
    "Parameter sweeps showing simulated efficiency vs radial confinement
    parameter, h_R").

All figures are saved as PNG files into ./outputs_plots/ (created if
missing) rather than shown interactively, so this runs headlessly.
"""

import os
import copy
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = "outputs_plots"


def _ensure_outdir():
    os.makedirs(OUTDIR, exist_ok=True)


def plot_transient(sol, title="", filename="transient.png"):
    _ensure_outdir()
    t = sol.t
    ne, ng, Te = sol.y

    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    axes[0].plot(t * 1e3, ne, color="tab:blue")
    axes[0].set_ylabel(r"$n_e$ [m$^{-3}$]")
    axes[0].set_yscale("log")

    axes[1].plot(t * 1e3, ng, color="tab:orange")
    axes[1].set_ylabel(r"$n_g$ [m$^{-3}$]")
    axes[1].set_yscale("log")

    axes[2].plot(t * 1e3, Te, color="tab:green")
    axes[2].set_ylabel(r"$T_e$ [eV]")
    axes[2].set_xlabel("time [ms]")

    fig.suptitle(title)
    fig.tight_layout()
    path = os.path.join(OUTDIR, filename)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  saved {path}")


def plot_hR_sweep(base_params, h_R_values=None, filename="hR_sweep.png"):
    """
    Sweep h_R (radial confinement parameter) and plot resulting efficiency,
    reproducing the structure of Wachs & Jorns Fig. 2a.
    """
    from models.jorns_magnetic_nozzle import rhs as jorns_rhs, performance as jorns_performance
    from solver import run_to_steady_state

    _ensure_outdir()
    if h_R_values is None:
        h_R_values = np.linspace(0.1, 0.9, 9)

    effs = []
    y0 = [1e16, 1e18, 3.0]   # warm-start seed; updated after each converged point
    for hR in h_R_values:
        p = copy.deepcopy(base_params)
        p.h_R = float(hR)
        try:
            _, y_ss, converged = run_to_steady_state(jorns_rhs, y0, p, t_span=(0.0, 2e-3), verbose=False)
            perf = jorns_performance(y_ss, p)
            effs.append(perf["efficiency"] * 100.0 if converged else np.nan)
            if converged:
                y0 = list(y_ss)   # warm-start next sweep point from this solution
        except Exception:
            effs.append(np.nan)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(h_R_values, effs, "o-", color="tab:red")
    ax.set_xlabel(r"Radial confinement parameter, $h_R$")
    ax.set_ylabel("Simulated thruster efficiency [%]")
    ax.set_title("Efficiency vs. radial confinement\n(cf. Wachs & Jorns, AIAA-2021-3382, Fig. 2a)")
    fig.tight_layout()
    path = os.path.join(OUTDIR, filename)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  saved {path}")
