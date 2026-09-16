# 0D Global Plasma Model for ECR Thrusters

A zero-dimensional (volume-averaged) global plasma discharge model for
Electron Cyclotron Resonance (ECR) electric-propulsion thrusters,
cross-referencing two distinct research lines:

1. **Jorns group (University of Michigan, PEPL)** — magnetic-nozzle ECR
   thruster architecture (no grids; ions accelerate/detach through an
   expanding magnetic nozzle):
   > Wachs, B. N. and Jorns, B. A., *"Optimization of an ECR Thruster using
   > Single, Two Frequency, and Pulsed Waveforms,"* AIAA Propulsion and
   > Energy Forum, **AIAA 2021-3382**, Aug. 2021.
   > (full text fetched from `pepl.engin.umich.edu/pdf/2021_AIAA_PE_Wachs.pdf`;
   > Eqs. (1)-(7) of that paper are implemented verbatim in
   > `models/jorns_magnetic_nozzle.py`.)

2. **Yamamoto group (Kyushu University / ISAS-JAXA μ-series lineage)** —
   gridded ECR ion thruster architecture (screen + accelerator grids
   electrostatically extract the ion beam; no magnetic nozzle):
   > Yamamoto, N., Kondo, S., Chikaoka, T., Nakashima, H., and Masui, H.,
   > *"Effects of magnetic field configuration on thrust performance in a
   > miniature microwave discharge ion thruster,"* **J. Appl. Phys. 102**,
   > 123304 (2007).
   > Dey, I., Toyoda, Y., Yamamoto, N., and Nakashima, H., *"Development of
   > a miniature microwave electron cyclotron resonance plasma ion thruster
   > for exospheric micro-propulsion,"* **Rev. Sci. Instrum. 86**, 123505
   > (2015).
   > and the broader μ1/μ10/MIPS lineage (Koizumi & Kuninaka 2010; Tani,
   > Yamashita, Tsukizaki, Nishiyama, Kuninaka, *Acta Astronautica* (2020);
   > Nishiyama et al., *Acta Astronautica* 166 (2020) — Hayabusa2 μ10 flight
   > operation).

Both architectures are built on the same underlying 0D global-model
machinery pioneered in Lieberman & Lichtenberg, *Principles of Plasma
Discharges and Materials Processing*, 2nd ed. (Wiley, 2005) — the textbook
both research lines cite as their common numerical ancestor — which is why
this code factors out a shared core (`models/base_model.py`,
`cross_sections.py`, `wall_losses.py`, `geometry.py`) and specializes only
the loss/extraction terms per architecture.

## IMPORTANT — what "research-accurate" means here

- The **Jorns/Wachs equations** (ion balance, neutral balance, power
  balance, h-factor bookkeeping, Bohm speed, sheath potential, Mach-number
  detachment) are transcribed **directly** from the fetched AIAA-2021-3382
  paper text and implemented as-is in `models/jorns_magnetic_nozzle.py`.
  This part is a faithful reproduction of a specific published equation set.
- The **Yamamoto-lineage gridded model** (`models/yamamoto_gridded.py`) is
  a **research-grounded synthesis**, not a verbatim transcription of one
  paper's internal equations (the exact closed-form ion/neutral balance
  used internally by ISAS's μ-series discharge-chamber papers was not
  available in full-text form during development — several were
  paywalled/robots-blocked). It applies the *same* standard global-model
  control-volume bookkeeping (Bohm flux, h-factors, sheath energy loss) that
  is standard across this entire literature (see e.g. Wang & Lafleur,
  *J. Appl. Phys.* 139, 083306 (2026), a contemporary gridded-ion-thruster
  global model using identical machinery), adapted for grid extraction
  instead of magnetic-nozzle detachment, and reports the same figures of
  merit the Yamamoto/Tsukizaki papers use (beam current, mass utilization,
  discharge loss in eV/ion).
- **Geometry and some h-factors/voltages for both example cases in
  `main.py` are illustrative**, not exact published thruster dimensions —
  the source papers either don't state them numerically in the fetched
  text, or (for ISAS hardware) are not fully public. Replace them with your
  own values for a publication-grade run. Anywhere a number is an
  assumption rather than a transcribed value, it is flagged in a comment.
- Ionization/excitation rate coefficients default to a **convenience
  analytic fit** for xenon (`cross_sections.K_iz_fit` /`K_ex_fit`). For a
  fully research-grade run, supply real tabulated cross sections (e.g.
  from [LXCat](https://fr.lxcat.net/)) via
  `cross_sections.load_lxcat_cross_section` +
  `RateCoefficientModel(species, iz_table=..., ex_table=...)`, which
  numerically integrates the Maxwellian average exactly as both research
  groups do internally.

## Project layout

```
ecr_global_model/
├── constants.py          physical constants, unit conversions
├── species.py            propellant properties (Xe/Kr/Ar)
├── cross_sections.py     Maxwellian-averaged rate coefficients (K_iz, K_ex)
├── geometry.py            chamber geometry (radius, length, grid open area)
├── wall_losses.py         Bohm velocity, sheath potential, h-factors
├── models/
│   ├── base_model.py               shared params dataclass + common physics
│   ├── jorns_magnetic_nozzle.py    Wachs & Jorns AIAA-2021-3382 model
│   └── yamamoto_gridded.py         Yamamoto-lineage gridded ion thruster model
├── solver.py              stiff transient integration + steady-state polish
├── plotting.py            transient + h_R parameter-sweep plots
├── main.py                 <-- RUN THIS. Runs both example cases end-to-end.
├── validation/
│   └── test_cases.py      internal consistency checks (mass/power balance closure)
├── requirements.txt
└── README.md              (this file)
```

## Installation & running

```bash
pip install -r requirements.txt
python3 main.py
```

This prints steady-state plasma properties (`n_e`, `n_g`, `T_e`) and
performance metrics (thrust, Isp, efficiency, mass utilization, and for the
gridded case, beam current and discharge loss) for both example cases, and
writes plots to `./outputs_plots/`:

- `jorns_transient.png` — time history of `n_e`, `n_g`, `T_e` approaching
  steady state (cf. Wachs & Jorns Fig. 1).
- `jorns_hR_sweep.png` — efficiency vs. radial confinement parameter `h_R`
  (cf. Wachs & Jorns Fig. 2a — efficiency should decrease monotonically as
  `h_R` increases, i.e. as radial confinement worsens).
- `yamamoto_transient.png` — same transient plot for the gridded case.

## Verifying it

Run the validation suite:

```bash
python3 validation/test_cases.py
```

This checks (not against experimental data, but internal self-consistency):
rate-coefficient positivity/monotonicity, steady-state particle-balance
closure (`dn_e/dt ≈ 0`, `dn_g/dt ≈ 0`), steady-state power-balance closure
(`P_abs = P_c + P_i + P_e [+ P_eT]` to <0.1% residual), and that efficiency
and mass-utilization are bounded in `[0, 1]` (catching, e.g., the exact bug
that was found and fixed during development — see below).

## Known modeling limitations (stated by the source papers themselves)

Directly from Wachs & Jorns, Sec. II: *"The model used here assumes fully
uniform plasma and neutral densities, Maxwellian electrons, a single ion
detachment point, and does not account for multiply charged ions, gas
heating, or multi-step ionization. As such, we use it for studying
performance trends and not absolute thrust predictions."* The same caveats
apply to the gridded model in this codebase, plus: no space-charge
(Child-Langmuir/perveance) limit on grid ion extraction is modeled (the
extracted current is Bohm-flux-limited only, which will overestimate beam
current at high density / low voltage relative to a real perveance-limited
grid).

## A note on a bug caught during development

An earlier version of `models/yamamoto_gridded.py::performance()` computed
wall-plug efficiency using only the microwave discharge power, omitting the
separate screen/accelerator grid beam-supply power (`I_beam * V_beam`).
Because gridded ion thrusters use two independent power supplies, and the
beam-supply power is often *larger* than the discharge power, this produced
an unphysical efficiency greater than 1. It was caught by exactly the kind
of check now automated in `validation/test_cases.py::test_yamamoto_performance_bounds`
and fixed by including `beam_power = I_beam * V_beam` in `total_power`. This
is left documented here deliberately, since "run it and see if the numbers
look sane" is a real and necessary part of using any global model like this
one.
