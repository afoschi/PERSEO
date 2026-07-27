# PERSEO - PERturbed Stellar Elements and Orbits.

Osculating orbital elements for stars orbiting a supermassive black hole under
a range of perturbing potentials.

This package integrates the Gauss planetary equations to compute the evolution
of the osculating elements (p, e, i, omega, Omega) over one orbit, for a test
particle subject to relativistic corrections, extended mass distributions, and
boson clouds. It was developed for orbital-dynamics studies of the S-stars at
the Galactic Center.

## Installation

Requires Python 3.9+, numpy and scipy.

    git clone https://github.com/[user]/[repo].git
    cd [repo]
    pip install -r requirements.txt

## Quick start

    import numpy as np
    from constants import G, c, Msun
    from relativistic import make_onePN
    from osculating_equations import integrate_orbit
    from postprocess import element_changes

    M = 4.29e6*Msun              # Sgr A* mass [kg]
    a = 1.5e14                   # semi-major axis [m]
    e = 0.88
    elements0 = (a*(1 - e**2), e, 2.0, 1.0, 3.0)   # p, e, i, omega, Omega

    perturbations = [make_onePN(M)]
    f, y = integrate_orbit(elements0, perturbations, M)

    d = element_changes(y, elements0, R0=8.28e3*pc_to_m)
    print("pericentre advance over one orbit:", d['omega'][-1], "deg")

## Available models

Each model is created by a factory that captures its physical parameters and
returns a callable `(f, elements) -> (R, S, W)`, the perturbing acceleration
components in the orbital frame. Models compose by superposition: pass a list
and the contributions are summed.

| Model | Factory | Module |
|---|---|---|
| 1PN Schwarzschild precession | `make_onePN` | `relativistic.py` |
| 1.5PN Lense-Thirring (spin) | `make_spin` | `relativistic.py` |
| Plummer sphere | `make_plummer` | `extended_mass.py` |
| Uniform ring | `make_ring` | `extended_mass.py` |
| Power-law surface density disk | `make_powerlaw_disk` | `extended_mass.py` |
| Scalar boson cloud | `make_scalar_cloud` | `boson_clouds.py` |
| Vector boson cloud | `make_vector_cloud` | `boson_clouds.py` |

Example of a combined model:

    perturbations = [
        make_onePN(M),
        make_spin(M, chi=1.0, i_spin=np.pi/2, beta_spin=0.0),
        make_powerlaw_disk(Mdisk, r_min, r_max, ibh=np.pi/2, beta=0.0, alpha=2.0),
    ]

## Conventions and units

All quantities are in SI units unless stated otherwise. The convention used for the angles is reported in Foschi et al. [Figure 1, arXiv] 

## Adding a new model

Write a factory returning a callable with signature `(f, elements)` that
returns the three acceleration components:

    def make_my_model(param1, param2):
        def my_model(f, elements):
            p, e, i, omega, Omega = elements
            ...
            return R, S, W
        return my_model

No changes elsewhere are needed.

## Tests

    python -m pytest tests/

## Citation

If you use this code, please cite [paper reference] and the Zenodo record:
[DOI badge once minted].

## License

MIT — see LICENSE.