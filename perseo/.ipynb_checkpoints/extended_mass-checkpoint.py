import numpy as np
from .constants import G
from scipy.special import ellipk, ellipe
from scipy.integrate import quad_vec



def power_law_density(radius, Mdisk, r_min, r_max, alpha):
    
    """ Computes the power law surface density and the mass of the disk between r_min and r """
    
    if (alpha == 2):
        B = Mdisk/(2*np.pi*np.log(r_max/r_min))
        Mr = 2*np.pi*B*np.log(radius/r_min) 
    else:
        B = Mdisk*(alpha-2)/(2*np.pi*(r_min**(2-alpha) - r_max**(2-alpha)))
        Mr = 2*np.pi*B*(r_min**(2-alpha) - radius**(2-alpha))/(alpha -2) 
    
    sigma = B/(radius)**(alpha)
                              
    return sigma, Mr


def elliptic(m):
    
    """ Complete elliptic integrals K(m), E(m).

    Uses the parameter convention m = k^2

    m : float or array
        Elliptic parameter. """
    
    K = ellipk(m)
    E = ellipe(m)
    
    return K, E


def make_plummer(r0, rho0):
    
    """ Perturbation from a Plummer density distribution [from Heißel et al. 2022]
        r0 : float
             Scale parameter [m] (usually r0 = 0.012 pc) 
        rho0: float
             Density at r = r0 [kg/m^3] (for r0 = 0.021 pc and rho0 = 1.69e-10 [kg/m^3], one recovers the upper limit within S2 orbit """
    
    def plummer(f, elements):
        p, e, i, omega, Omega = elements
        r = p/(1 + e*np.cos(f))
    
        R = -4*np.pi*G*rho0*r*r0**3/(3*(r**2 + r0**2)**(3/2))
        W = np.zeros_like(R)
        S = np.zeros_like(R)
        return R, S, W
    return plummer

def make_ring_from_sigma(sigma, R_ring, ibh, beta):
    
    """ Uniform ring specified by linear surface density rather than total mass with sigma = Mring/(2*pi*R_ring). """
    
    def uniform_ring(f, elements):
        p, e, i, omega, Omega = elements
        gamma = f - beta + omega 
        r = p/(1 + e*np.cos(f))
        
        sg, cg = np.sin(gamma), np.cos(gamma)
        si, ci = np.sin(ibh), np.cos(ibh)
        if si <= 0.01:
            si = 0.0
  
        d = r*np.sqrt(ci**2*cg**2 + sg**2)
        zBH = r*cg*si

        pa = np.sqrt((d + R_ring)**2 + zBH**2)
        qa = np.sqrt((d - R_ring)**2 + zBH**2)
        ka = 4*R_ring*d/pa**2

        K, E = elliptic(ka)
    
        #random choice that does not influence the results 
        OBH = 0.0

        Ap = 2*G*sigma*R_ring/pa*(K/d**2 +(E/qa**2)*(1 - (R_ring**2 + zBH**2)/d**2))
        Aq = 4*G*sigma*R_ring*E/(pa*qa**2)
    
        xBH = -r*(ci*cg*np.sin(OBH) + np.cos(OBH)*sg)
        yBH = r*(ci*cg*np.cos(OBH) - np.sin(OBH)*sg)
    
        Ax, Ay, Az = -Ap*xBH, -Ap*yBH, -Aq*zBH
    
        Rcomp = -sg*(Ax*np.cos(OBH) + Ay*np.sin(OBH)) + cg*(Az*si + ci*(Ay*np.cos(OBH) - Ax*np.sin(OBH)))
        W = Az*ci + si*(-Ay*np.cos(OBH) + Ax*np.sin(OBH))
        S = -(Az*si + ci*(Ay*np.cos(OBH) - Ax*np.sin(OBH)))*sg - (Ax*np.cos(OBH) + Ay*np.sin(OBH))*cg
        return Rcomp, S, W
    return uniform_ring


def make_ring(Mring, R_ring, ibh, beta):
     
    """ Perturbation due to a uniform ring 
        Mring : float
             Mass of the ring [kg]  
        R: float
            Radius of the ring [m]
        ibh: float
            Inclination between the ring z axis and the orbital axis z_orb [rad]
        beta: float
            Angle between the orbital line of nodes and the projection of the ring z axis onto the orbital plane [rad] """
    
    return make_ring_from_sigma(Mring/(2*np.pi*R_ring), R_ring, ibh, beta)


def make_powerlaw_disk(Mdisk, r_min, r_max, ibh, beta, alpha, log_spacing=True):
    
    """ Perturbation due to a disk with power law surface density Sigma = B * r^{-alpha}
        Mdisk : float
             Mass of the disk between r_min and r_max [kg]  
        r_min: float
            Inner radius of the disk [m]
        r_max: float
            Outer radius of the disk [m]
        ibh: float
            Inclination between the disk z axis and the orbital axis z_orb [rad]
        beta: float
            Angle between the orbital line of nodes and the projection of the disk z axis onto the orbital plane [rad] 
        alpha: float
            Power law slope """
    
    
    def power_law_disk(f, elements):
        p, e, i, omega, Omega = elements
    
        def components(radius):
            sigma, _ = power_law_density(radius, Mdisk, r_min, r_max, alpha)
            ring = make_ring_from_sigma(sigma, radius, ibh, beta)
            return ring(f, elements)

        if log_spacing:
            lo, hi = np.log(r_min), np.log(r_max)
            integrand = lambda u: np.asarray(components(np.exp(u)))*np.exp(u)
        else:
            lo, hi = r_min, r_max
            integrand = lambda x: np.asarray(components(x))

        result, _ = quad_vec(integrand, lo, hi, epsrel=1e-8)
        return tuple(result)
    return power_law_disk


