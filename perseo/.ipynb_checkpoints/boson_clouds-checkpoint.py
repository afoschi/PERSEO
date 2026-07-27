import numpy as np
from .constants import G


def make_scalar_cloud(M, Mcloud, alpha):
    
    """ Perturbation induced by a scalar cloud placed on the equatorial plane of the BH (theta = pi/2).  
        M : float
            Mass of the central object [kg] 
        alpha : float
            Coupling between gravity and scalar field [dimensionless] """
    mu = G*Mcloud/M
    
    def scalar_cloud(f, elements):
        p, e, i, omega, Omega = elements
        r = p/(1 + e*np.cos(f))
        p1, p2, dp1, dp2 = perturbative_functions(r, alpha, M)
    
        R = mu*(np.sin(i)**2*np.sin(omega +f)**2*dp2 + dp1)
        S = mu*(np.sin(i)**2*(e*np.cos(f)+1)*np.sin(2*(f + omega))*p2)/p
        W = mu*(np.sin(2*i)*(e*np.cos(f) + 1)*np.sin(f + omega)*p2)/p 
        return R, S, W
    return scalar_cloud

def make_vector_cloud(M, Mcloud, alpha):
    
    """ Perturbation induced by a vector cloud.  
        M : float
            Mass of the central object [kg] 
        Mcloud : float
            Mass of the cloud [kg]
        alpha : float
            Coupling between gravity and vector field [dimensionless] """

    mu = G*Mcloud/M

    def vector_cloud(f, elements):
        p, e, i, omega, Omega = elements
        r = p/(1 + e*np.cos(f))

        expfac = np.exp(-2*r*alpha**2/M)
        
        Rcomp = - mu*M/r**2 + mu*expfac/(M*r**2) * (M**2 + 2*M*r*alpha**2 + 2*r**2*alpha**4)
        S = np.zeros_like(Rcomp)
        W = np.zeros_like(Rcomp)
        return Rcomp, S, W
    return vector_cloud


def perturbative_functions(r, alpha, M):
    alpha_2 = alpha**2
    alpha_4 = alpha**4
    alpha_6 = alpha**6
    
    m = (G*M)/r
    
    exp = np.exp(-alpha**2/m)
    
   
    P1 = m + (3*m**3)/(alpha_4) - exp*((3*m**3)/alpha_4 + 3*m**2/alpha_2 + 5*m/2 + 5*alpha_2/4 \
                                      + (3*alpha_4)/(8*m) + alpha_6/(16*m**2))
         
        
    P2 = -9*(m**3)/(alpha_4) + exp*((9*m)/2 + (9*m**3)/(alpha_4) + (9*m**2)/(alpha_2)\
                                      + (3*alpha**2)/2 + (3*alpha_4)/(8*m) \
                                      + (alpha_6)/(16*m**2))
    
    
    dP1 = - m/r - 9*m**3/(alpha_4*r) + (exp/(16*r))*(88*m + 144*m**3/alpha_4 + 144*m**2/alpha_2 \
                                                    + 40*alpha_2 + 14*alpha_4/m + 4*alpha_6/m**2 \
                                                    + alpha_6*alpha_2/m**3)
        
    dP2 = (27*m**3)/(alpha_4*r) - (exp/(16*r))*(216*m + 432*m**3/alpha_4 + 432*m**2/alpha_2 \
                                               + 72*alpha_2 + 18*alpha_4/m + 4*alpha_6/m**2 \
                                               + alpha_6*alpha_2/m**3)
    
    return  P1, P2, dP1, dP2