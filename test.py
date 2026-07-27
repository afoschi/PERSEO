"""
Basic validation tests for PERSEO perturbing accelerations.

Cross-checks against:
  - Plummer enclosed mass: GRAVITY (2024) & Heissel et al. (2022) upper limit within the S2 orbit (1000 Msun < M < 4000 Msun).
  - 1PN Schwarzschild periapsis advance: GRAVITY (2020) / Heissel et al. (2022), ~12 arcmin per orbit for S2.

Both make_plummer and make_onePN return perturbing accelerations (R, S, W) in the
RSW frame as functions of true anomaly f, so the tests reconstruct the physical
quantities (enclosed mass, per-orbit precession) from those components.

Run with:  python -m pytest tests/ -v
"""

import numpy as np
import pytest

from perseo.extended_mass.make_plummer import make_plummer
from perseo.relativistic.make_onePN import make_onePN
from perseo.osculating_equations import integrate_orbit, element_changes

# --- Physical constants (SI) -------------------------------------------------
G = 6.674e-11          # m^3 kg^-1 s^-2
c = 2.99792458e8       # m s^-1
Msun = 1.98892e30      # kg
pc = 3.0857e16         # m
au = 1.495978707e11    # m
arcsec = np.pi / (180.0 * 3600.0)   # rad
deg = 0.0174533 # from deg to rad

# --- Sgr A* parameters --------------------------------------------------
M_BH = 4.297e6 * Msun

# S2 orbit [from GRAVITY (2024)]

a_S2 = 1.5472e14      # m 
e_S2 = 0.8844
i_S2 = 134.7*deg
omega_S2 = 66.3*deg + np.pi
Omega_S2 = 228.2*deg + np.pi #note the shift by a factor pi as explained in Foschi et al. (2026)
p_S2 = a_S2 * (1 - e_S2**2)

# -----------------------------------------------------------------------------
# Plummer: enclosed-mass cross-check
# -----------------------------------------------------------------------------
def plummer_enclosed_mass(r, r0, rho0):
    """ Analytic enclosed mass of a Plummer sphere at radius r.

    M(<r) = (4/3) pi rho0 r0^3 * r^3 / (r^2 + r0^2)^(3/2)

    This is the closed-form integral of the Plummer density whose radial
    acceleration make_plummer returns as R = -G M(<r) / r^2. """
    
    return (4.0 / 3.0) * np.pi * rho0 * r0**3 * r**3 / (r**2 + r0**2)**1.5


def test_plummer_R_matches_enclosed_mass():
    """ The radial acceleration from make_plummer must equal -G M(<r)/r^2
    for the analytic Plummer enclosed mass, at an arbitrary point on the orbit. """
    
    r0 = 0.021 * pc
    rho0 = 1.69e-10  # kg/m^3  (docstring config: S2 upper limit)

    plummer = make_plummer(r0, rho0)

    # Build elements with a semi-latus rectum so that at f=0, r = p/(1+e).

    elements = (p_S2, e_S2, 0.0, 0.0, 0.0)

    f = 0.0
    R, S, W = plummer(f, elements)

    r = p_S2 / (1 + e_S2 * np.cos(f))
    R_expected = -G * plummer_enclosed_mass(r, r0, rho0) / r**2

    assert np.isclose(R, R_expected, rtol=1e-10)
    assert S == 0.0 and W == 0.0


def test_plummer_enclosed_mass_within_S2_is_upper_limit():
    """ With the docstring config (r0=0.021 pc, rho0=1.69e-10), the enclosed
    Plummer mass within the S2 orbit should sit at the ~ (1-5) 10^3 Msun upper
    limit reported for extended mass inside S2 [GRAVITY (2024)]. """
    
    r0 = 0.021 * pc
    rho0 = 1.69e-10

    r_apo = a_S2 * (1 + e_S2)   # apoapsis radius
    M_enc = plummer_enclosed_mass(r_apo, r0, rho0) / Msun

    assert 1e3 < M_enc < 5e3


# -----------------------------------------------------------------------------
# 1PN Schwarzschild: periapsis advance cross-check
# -----------------------------------------------------------------------------
#
# make_onePN returns (R, S, W) as functions of true anomaly. The periapsis
# advance per orbit comes from integrating Gauss's planetary equation for omega
# over one revolution. The analytic 1PN result must reproduce:
#
#     Delta_omega = 6 pi G M / (c^2 a (1 - e^2)) [rad per orbit]
#                 ~ 12.1 arcmin/orbit [for S2] [GRAVITY (2020)]
#


def onePN_precession_analytic(M_BH, p_S2):
    """ Analytical-form 1PN periapsis advance per orbit [rad]. """
    return 6.0 * np.pi * G * M_BH / (c**2 * p_S2)


def test_onePN_integrated_advance_matches_analytic():
    """ Integrate make_onePN(M) over one S2 orbit and compare the resulting
    periapsis advance to the analytic 6 pi GM / (c^2 p) value. """

    perturbations = [make_onePN(M_BH)]

    elements0 = (p_S2, e_S2, i_S2, omega_S2, Omega_S2)
    f, y = integrate_orbit(elements0, perturbations, M_BH)

    dw = element_changes(f, y, elements0, angle_units='rad')[3]
    assert np.isclose(dw, onePN_precession_analytic(M_BH, p_S2),
                          rtol=1e-2)