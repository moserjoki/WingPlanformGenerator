import math
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from sizing import *

M_cruise = 0.85 # []
m_MTOW = 114960 # [kg] Maximum Take Off Weight from Class I Weight Estimation . 
AR = 6.9 # [] aspect ratio
g = 9.81 # [m/s^2] gravitational constant
ρ_cruise = 0.44 # [kg/m^3] density at cruise altitude 
V_cr = 253 # [m/s] Cruise spseed 

# Values from chosen airfoil SC(2)-0414
c_l_alpha = 6.6954 # []
c_d0 = 0.00551 # [] 
tau = 0.48 # []
c_l_max = 1.528 # [] 

# Create objects used for subsequent calculations
airfoil = Airfoil("airfoils/NASA SC(2)-0414.dat", c_l_alpha, c_d0, tau, c_l_max) 
wing = WingSizing(m_MTOW, AR, V_cr, M_cruise, ρ_cruise, airfoil)
cruise_matching_diagram = MatchingDiagram(m_MTOW, AR, V_cr, M_cruise, ρ_cruise)


# Dummy values for first matching diagram calculation
C_L_max_take_off_cur = 2.1 # [] Cl during take-off
C_L_max_landing_cur = 2.4 # [] Cl during landing



# Iterating on the design based on the initial given configuration until it converges
for i in range(5):
    print(f"Iteration {i+1}")
    S_w_cur  = cruise_matching_diagram.compute(C_L_max_take_off_cur, C_L_max_landing_cur)
    wing.planform_sizing(S_w_cur)
    wing.aileron_sizing(C_L_max_landing_cur)

    # DATCOM methode 
    C_L_C_l_max_ratio = 0.95 # [] From graph
    delta_C_L_max = 0 # [] From graph
    C_L_max_clean = wing.DATCOM_C_L_max_clean(airfoil, C_L_C_l_max_ratio, delta_C_L_max)

    C_L_max_take_off_cur, C_L_max_landing_cur = wing.HLD_sizing(C_L_max_clean)
 
cruise_matching_diagram.plot()


#FUSELAGE PARAMETERS:
l_fus=52.43 #RANDOM INITIAL VALUE
X_cg_aft=21.98 #RANDOM INITIAL VALUE

wing.empenage_sizing(l_fus, X_cg_aft, True)

# Wing fuel tank parameters
xc_f1 = 0.2 # [] cord ratio from which onwards fuel is stored in the cross section 
xc_f2 = 0.65 # [] cord ratio until which fuel is stored in the cross section
dxc = 0.001 # [] iteration steps for cross section integration
b_f1 = (wing.b/2)*0 # [m] span start position from which onwards the fuel is stored in the wing
b_l1 = (wing.b/2)*0.3 # [m] span start position of landing gear omission region where no fuel can be stored
b_l2 = (wing.b/2)*0.3 + 2 # [m] span end position of landing gear where no fuel can be stored 
b_f2 = (wing.b/2)*0.65 # [m] span end position until which the fuel is stored in the wing 
dbf = 0.01 # [m] iteration steps for span section integration

# Wing fuel tank volume determination
S_ref_chord = airfoil.cross_section(xc_f1, xc_f2, dxc, True)

V_available_volume = wing.volume_section(b_f1, b_f2, S_ref_chord, dbf)
V_landing_gear = wing.volume_section(b_l1, b_l2, S_ref_chord, dbf) 

V_fuel_half_wing = V_available_volume - V_landing_gear 
V_fuel_total_wing = 2*V_fuel_half_wing

print(f"Cross Section: {S_ref_chord}")
print(f"Volume: {V_fuel_total_wing} m^3")
wing.plot()