# Conventions: 
# Main program has variables that change throughout the iterations
# Parameters has variables that are either pre-set or have to be changed manually by updating information from graphs based on design decisions

# Cruise Requirement 
M_cruise = 0.85 # []
ρ_cruise = 0.44 # [kg/m^3] density at cruise altitude 
V_cruise = 253 # [m/s] Cruise spseed 

# Atmospheric constants
ρ_sea_level = 1.22522568 # [kg/m^3] Density at sea level

# DATCOM methode 
datcom_C_L_C_l_max_ratio = 0.95 # [] From graph
datcom_delta_C_L_max = 0 # [] From graph

# Values from chosen airfoil SC(2)-0414
airf_c_l_alpha = 6.6954 # []
airf_c_d0 = 0.00551 # [] 
airf_tau = 0.48 # []
airf_c_l_max = 1.528 # [] 

# General constants
g = 9.81 # [m/s^2] gravitational acceleration

# Fuselage Parameters
d_fuselage = 5.2 # [m] Fuselage diameter 
l_fus = 44.427 # [m] Fuselage length


# Mach Drag Divergence
mdd_k_a = 0.935             # [] technology factor for all super critical airfoils
mdd_t_c_streamwise = 0.14   # [] t/c 

# Aileron Parameters
ail_b2_percent = 0.75 # [] Most outward position of ailerons (not to tip to minimize alerion reversal)
ail_delta_a_up = 17 # [°] Deflection angle of the upward moving aileron
ail_f_differential_aileron = 0.75 # [] Ratio of the deflection angle of the upward moving to the downward moving aileron 
ail_delta_bank_req = 60  # [°] Roll requirement based on CS25 
ail_delta_t_req = 7 # [s] Max time required to fulfill the roll req based on CS25

# High Lift Device Parameters
hld_S_wfLE_to_S = 0.8 # [] Leading Edge Fraction Flapped Area:
hld_delta_C_l_LE = 0.4 # !G [] Slat
hld_delta_C_l_TE_take_off = 1.88 # !G [] Single slotted Fowler flap 
hld_delta_C_l_TE_landing = 1.94 # !G [] Single slotted Fowler flap 

# Empenage Parameters
# Selected components. All values were selected to be in the middle of acceptable range:
empg_Vh = 1.01
empg_Vv = 0.079
empg_AR_v = 1.5
empg_AR_h = 4
empg_taper_v = 0.5
empg_taper_h = 0.65

# Wing fuel tank parameters
fuel_xc_1 = 0.2 # [] cord ratio from which onwards fuel is stored in the cross section 
fuel_xc_2 = 0.65 # [] cord ratio until which fuel is stored in the cross section
fuel_dxc = 0.001 # [] iteration steps for cross section integration
fuel_b_f1 = 0 # [m] span start position from which onwards the fuel is stored in the wing
fuel_b_l1 = 0.3 # [m] span start position of landing gear omission region where no fuel can be stored
fuel_delta_l1_l2 = 2 # [m] delta span from span start where no fuel can be stored because of landing gear
fuel_b_f2 = 0.65 # [m] span end position until which the fuel is stored in the wing 
fuel_dbf = 0.01 # [m] iteration steps for span section integration