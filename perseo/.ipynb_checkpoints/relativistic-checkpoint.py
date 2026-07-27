import numpy as np
from .constants import G, c


def make_onePN(M):
    """ First Post Newtonian correction to Newtonian 2-body motion.
        M : float
            Mass of the central object [kg] """
    
    def onePN(f, elements):
        p, e, i, omega, Omega = elements
        R = (G**2*M**2/(c**2*p**3))*(1+e*np.cos(f))**2*(3*(e**2+1)+2*e*np.cos(f) - 4*e**2*np.cos(f)**2)
        S = (4*G**2*M**2/(c**2*p**3))*(1 + e*np.cos(f))**3*e*np.sin(f)
        W = np.zeros_like(R)
        return R, S, W
    return onePN


def make_spin(M, chi, i_spin, beta_spin):
    """ 1.5 Post Newtonian correction to Newtonian 2-body motion.
        M : float
            Mass of the central object [kg]
        chi : float
            Magnitude of the spin (between 0 and 1)
        i_spin : float
            Inclination between BH spin axis and z-axis of the orbit z_orb [rad]
        beta_spin : float
            Angle between the orbital line of nodes and the projection of the BH spin axis onto the orbital plane [rad] """
    
    def spin(f, elements):
        p, e, i, omega, Omega = elements
        eps = (G*M/(c**2*p))**(3/2)
        phi = omega + f
    
        R = 2*eps*G*M*chi*(1+e*np.cos(f))**4 * np.cos(i_spin)/p**2
    
        S = -2*eps*G*M*chi*e*np.sin(f)*(1 + e*np.cos(f))**3*np.cos(i_spin)/p**2
        W = 2*eps*G*M*chi*(1 + e*np.cos(f))**3*np.sin(i_spin)/p**2 * (2*(1 + e*np.cos(f))*np.cos(beta_spin - phi) \
                                                                      + e*np.sin(f)*np.sin(beta_spin - phi))
        return R, S, W
    return spin