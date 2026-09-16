"""
geometry.py
-----------
Discharge-chamber geometry for a cylindrical ECR thruster source, common to
both reference architectures:

  - Jorns-group magnetic-nozzle ECR thruster (Wachs & Jorns, AIAA-2021-3382):
    cylindrical source of radius R, length L, feeding an open magnetic
    nozzle at the exit plane (area pi*R^2) -- no grids.

  - Yamamoto-group gridded ECR ion thruster (mu1/mu10/MIPS lineage,
    Yamamoto et al., J. Appl. Phys. 102, 123304 (2007)): cylindrical
    source of radius R, length L, with a screen-grid aperture of open
    area A_screen << pi*R^2 at the exit plane through which ions are
    extracted, and the remainder of the exit plane (back wall / screen
    grid solid area) acting as an additional wall-loss surface.

All areas below follow the same "control volume" bookkeeping used in
Wachs & Jorns Eq. (1)-(2): radial wall area 2*pi*R*L, and end-plate
("back wall" + "exit/screen") areas each pi*R^2.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class ThrusterGeometry:
    R: float                 # chamber radius [m]
    L: float                 # chamber length [m]
    screen_open_area_frac: float = 1.0
    # For the gridded (Yamamoto-style) architecture, the fraction of the
    # exit-plane area that is physically open (screen grid transparency).
    # Set to 1.0 for the magnetic-nozzle (Jorns-style) architecture, where
    # the entire exit plane is an open magnetic nozzle throat.

    @property
    def volume(self):
        """Chamber volume V = pi R^2 L [m^3]."""
        return np.pi * self.R ** 2 * self.L

    @property
    def area_radial(self):
        """Lateral (radial) wall area = 2 pi R L [m^2]."""
        return 2.0 * np.pi * self.R * self.L

    @property
    def area_end(self):
        """One end-plate area = pi R^2 [m^2] (back wall OR exit plane)."""
        return np.pi * self.R ** 2

    @property
    def area_screen_open(self):
        """Open (grid aperture / nozzle throat) area seen by escaping ions."""
        return self.area_end * self.screen_open_area_frac

    @property
    def area_screen_solid(self):
        """Solid (non-transparent) fraction of the exit end-plate, which
        acts as an additional ion/electron wall-loss surface (relevant for
        the gridded architecture; is zero when screen_open_area_frac = 1)."""
        return self.area_end * (1.0 - self.screen_open_area_frac)
