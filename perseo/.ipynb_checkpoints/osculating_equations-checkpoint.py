import numpy as np
from scipy.integrate import solve_ivp
from .constants import G, M, c, as_to_m

def gauss_rhs(f, state, perturbations, M):
    p, e, i, omega, Omega = state[:5]

    R = S = W = 0.0
    for pert in perturbations:
        dR, dS, dW = pert(f, state[:5])
        R += dR; S += dS; W += dW

    u = 1.0 + e*np.cos(f)
    n = np.sqrt(p/(G*M))

    dOmegadt = n*np.sin(omega + f)/u*(W/np.sin(i))
    domegadt = (1/e)*n*(-np.cos(f)*R + (2 + e*np.cos(f))*np.sin(f)*S/u
                        - e*(1/np.tan(i))*np.sin(omega + f)*W/u)
    didt = n*np.cos(f + omega)*W/u
    dpdt = 2*np.sqrt(p**3/(G*M))*S/u
    dedt = n*(np.sin(f)*R + (2*np.cos(f) + e*(1 + np.cos(f)**2))*S/u)

    # derived angles, defined through the element rates
    dvarpidt = domegadt + np.cos(i)*dOmegadt
    dThetadt = np.sin(omega)*didt - np.sin(i)*np.cos(omega)*dOmegadt
    dXidt    = -np.cos(omega)*didt - np.sin(i)*np.sin(omega)*dOmegadt

    dfdt = np.sqrt(G*M/p**3)*u**2 - dvarpidt
    dtdf = 1.0/dfdt

    rates = np.array([dpdt, dedt, didt, domegadt, dOmegadt,
                      dvarpidt, dThetadt, dXidt])
    return rates*dtdf



def integrate_orbit(elements0, perturbations, M, n_points=100, **kwargs):
    if not perturbations:
            raise ValueError("No perturbations applied; unperturbed orbit has constant elements.")
    state0 = np.concatenate([np.asarray(elements0, dtype=float), np.zeros(3)])

    f_span = (1e-10, 2*np.pi)
    f_eval = np.linspace(*f_span, n_points)
    sol = solve_ivp(gauss_rhs, f_span, state0, method = 'Radau',
                    t_eval = f_eval, args = (perturbations, M),
                    dense_output = True, rtol = 1e-6, atol = 1e-6, **kwargs)
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol.t, sol.y

def element_changes(f, y, elements0, angle_units='deg'):

    """ Change in osculating elements relative to their initial values.

    f : array
        True anomaly grid [rad].
    y : array, shape (5, N)
        Elements (p, e, i, omega, Omega) along the orbit.
    elements0 : sequence
        Initial elements, as passed to integrate_orbit (varpi0, Theta0, Xi0 = 0)
    as_to_m is defined as as_to_m = as_to_rad*R0 [m] with R0 = 2.55e20 m. 
   
    Returns
    -------
    dict with keys 'p', 'e', 'i', 'omega', 'Omega', 'varpi', 'Theta', 'Xi'  """
    
    p0, e0, i0, omega0, Omega0 = elements0

    
    scale = {'rad': 1.0,
             'deg': 180.0/np.pi,
             'arcsec': 180.0/np.pi*3600.0}[angle_units]
    
    return {
        'p':     (y[0] - p0)/as_to_m,
        'e':      y[1] - e0,
        'i':     (y[2] - i0)*scale,
        'omega': (y[3] - omega0)*scale,
        'Omega': (y[4] - Omega0)*scale,
        'varpi':  y[5]*scale,      # already a change, no subtraction
        'Theta':  y[6]*scale,
        'Xi':     y[7]*scale,
    }
    
