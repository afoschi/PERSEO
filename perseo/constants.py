G    = 6.67430e-11
c    = 299792458
Msun = 1.98841e30

years_to_s = 365.25*24*3600

s_to_kg   = c**3/G
m_to_kg   = c**2/G  
as_to_rad = 4.848136811095360e-06
degree_to_rad = 0.0174532925199433
pc_to_m = 3.08567758149e16 
m_to_au = 6.68e-12

conv2 = 1e6*Msun

as_to_kg = as_to_rad*(8277.09055007*pc_to_m)*m_to_kg
as_to_m = as_to_rad*(8277.09055007*pc_to_m)

e = 0.884429099282                           # Eccentricity


R0 = 8277.09055007*pc_to_m                   # Distance from the Galactic Center [m]
R0_kg = R0*m_to_kg                           # Distance from the Galactic Center [kg]
R0_new = R0_kg/conv2


M = 4.29701742727e6*Msun
M_new = M/conv2

tp = 2018.3789878
tp_s = 2018.3789878*years_to_s     # Periastron passage [s]
tp_kg = tp_s*s_to_kg               # Periastron passage [kg]
tp_new = tp_kg/conv2


a = 0.1249527719
a_rad = a*as_to_rad     # Semi-major axis [rad]
a_m = a_rad*R0                     # Semi-major axis [m] 
a_kg = a_m*m_to_kg                 # Semi-major axis [kg]
a_new = a_kg/M                     # Semi-major axis normalized by BH mass (adimensional, G.U.)


OmegaORB = 228.191510132*degree_to_rad      # Longitude of the periastron [rad]
omegaORB = 66.2689390128*degree_to_rad      # Argument of the periastron [rad]
iORB = 134.700204975*degree_to_rad          # Inclination [rad]



RASTAconv = years_to_s*s_to_kg/conv2



# priors on offsets 

y0_prior = -0.00057 
sigma_y0 = 0.00015

x0_prior = -0.000055
sigma_x0 = 0.00025

vz0_prior = 0
sigma_vz0 = 5

vx0_prior = 0.000032
sigma_vx0 = 0.000019

vy0_prior = 0.000063
sigma_vy0 = 6.6e-6

