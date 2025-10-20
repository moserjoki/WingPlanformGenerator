import math
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from sizing import *
from drag import *

m_MTOW = 114960 # [kg] Maximum Take Off Weight from Class I Weight Estimation . 
AR = 6.9 # [] aspect ratio

# Create objects used for subsequent calculations
airfoil = Airfoil("airfoils/NASA SC(2)-0414.dat") 
wing = WingSizing(m_MTOW, AR)
cruise_matching_diagram = MatchingDiagram(m_MTOW, AR)

# Dummy values for first matching diagram calculation
C_L_max_take_off_cur = 2.1 # [] Cl during take-off
C_L_max_landing_cur = 2.4 # [] Cl during landing
S_w_cur = 0

# Iterating on the design based on the initial given configuration until it converges
for i in range(5):
    print(f"Iteration {i+1}")
    S_w_cur  = cruise_matching_diagram.compute(C_L_max_take_off_cur, C_L_max_landing_cur)
    wing.planform_sizing(S_w_cur)
    wing.aileron_sizing(C_L_max_landing_cur)

    C_L_max_clean = wing.DATCOM_C_L_max_clean()
    C_L_max_take_off_cur, C_L_max_landing_cur = wing.HLD_sizing(C_L_max_clean)
 
#cruise_matching_diagram.plot()

X_cg_aft = 21.98 #RANDOM INITIAL VALUE
b_v, c_r_v, c_t_v, MAC_v, b_h, c_r_h, c_t_h, MAC_h = wing.empenage_sizing(X_cg_aft, True)

#wing.fuel_volume(airfoil)

#wing.plot()

# === INPUT PARAMETERS FOR EACH COMPONENT ===
# Flight conditions (same for all components)
flight_condition = 'cruise'  # 'takeoff' or 'cruise'

# Get atmospheric properties based on flight condition
atm_properties = get_atmospheric_properties(flight_condition)
rho = atm_properties['rho']
altitude = atm_properties['altitude']

# Flight conditions configuration:
if flight_condition == 'landing':
    V = prm.V_approach
    M = prm.M_approach
    flow_regime = 'subsonic'
    gear_deployed = True
    flaps_deflected = True
    flap_deflection_deg = 20.0  # Takeoff flap setting
else:  # cruise
    V = prm.V_cruise
    M = prm.M_cruise
    flow_regime = 'transonic'
    gear_deployed = False
    flaps_deflected = False
    flap_deflection_deg = 0.0   # No flap deflection in cruise

# Landing gear parameters
CD_gear = 0.015         # Constant drag coefficient for landing gear

# Flap parameters (slotted flaps)
c_f = 0.8               # Chord length of flap [m]
c = 2.0                 # Wing chord length [m]
S_flap = 10.0           # Wing area affected by flap [m²]

# Excrescence and leakage parameters
excrescence_percentage = 0.035  # 3.5% of total CD0 (2-5% range)

# Calculate wetted areas
S_wet_wing = Swet_wing(wing.b,wing.c_tip,wing.chord(d_fuselage/2),d_fuselage)
S_wet_VT = Swet_VT(b_v, c_r_v, c_t_v)
S_wet_HT = Swet_HT(b_h, c_r_h, c_t_h)
S_wet_fus = Swet_fus()

components = {
    'wing': {
        'l': wing.MAC, # [m] Wing MAC 
        'laminar_fraction': 0.1, # 10% laminar flow on wing
        'type': 'lifting_surface',
        't_c': prm.airf_t_c,
        'x_c_max': prm.airf_x_c_max, 
        'Lambda_max_deg': 29.69524,  # Sweep at max thickness [deg]
        'IF': 1.3,  # Wing-fuselage interference
        'wetted_area': S_wet_wing
    },
    'fuselage': {
        'l': prm.l_fus,          # Fuselage length [m]
        'laminar_fraction': 0.05,  # 5% laminar flow on fuselage
        'type': 'fuselage',
        'A_max': prm.A_max_fus, # Maximum cross-sectional area [m²]
        'IF': 1.0,  # Fuselage base drag, protuberances
        'wetted_area': S_wet_fus
    },
    'HT': {
        'l': MAC_h,           # [m] Horizontal Tail MAC 
        'laminar_fraction': 0.10,  # 5% laminar flow on HT
        'type': 'lifting_surface', 
        't_c': 0.10,        # Thickness to chord ratio
        'x_c_max': 0.3,     # Location of max thickness
        'Lambda_max_deg': 33.78,  # Sweep at max thickness [deg]
        'IF': 1.04,  # Tail-fuselage interference
        'wetted_area': S_wet_HT
    },
    'VT': {
        'l': MAC_v,           # [m] Vertical Tail MAC 
        'laminar_fraction': 0.10,  # 5% laminar flow on VT
        'type': 'lifting_surface',
        't_c': 0.10, # Thickness to chord ratio
        'x_c_max': 0.3, # Location of max thickness
        'Lambda_max_deg': 26.98,  # Sweep at max thickness [deg]
        'IF': 1.04,  # Tail-fuselage interference
        'wetted_area': S_wet_VT
    }
}

print("\n\n\n")
print("\n" + "=" * 120)
# Print flight condition information
print(f"FLIGHT CONDITION: {flight_condition.upper()}")
print(f"Altitude: {altitude}")
print(f"Air Density: {rho:.3f} kg/m³")
print(f"Velocity: {V:.1f} m/s")
print(f"Mach Number: {M:.2f}")
print(f"Flow Regime: {flow_regime}")
print("=" * 120)

# Calculate skin friction for each component
print("Component Drag Calculations:")
print("=" * 120)
print(f"{'Component':<12} {'Swet [m²]':<10} {'Cf':<8} {'FF':<8} {'IF':<8} {'Cf×FF×IF':<12} {'Cf×FF×IF×Swet':<15}")
print("-" * 120)

results = {}
total_aircraft_drag_component = 0

for component, params in components.items():
    # Calculate skin friction
    result = calculate_component_skin_friction(
        rho, V, params['l'], k, M, flow_regime, params['laminar_fraction'], flight_condition
    )
    
    # Calculate form factor based on component type
    if params['type'] == 'lifting_surface':
        FF = form_factor_lifting_surfaces(
            params['t_c'], params['x_c_max'], M, params['Lambda_max_deg']
        )
    elif params['type'] == 'fuselage':
        FF = form_factor_fuselage(params['l'], params['A_max'], result['Re_effective'])
    else:
        FF = 1.0  # Default form factor
    
    # Get interference factor and wetted area
    IF = params['IF']
    Swet = params['wetted_area']
    
    # Calculate final drag coefficient components
    Cf_FF_IF = result['Cf_total'] * FF * IF
    Cf_FF_IF_Swet = Cf_FF_IF * Swet
    
    # Accumulate total
    total_aircraft_drag_component += Cf_FF_IF_Swet
    
    result['form_factor'] = FF
    result['IF'] = IF
    result['wetted_area'] = Swet
    result['Cf_FF_IF'] = Cf_FF_IF
    result['Cf_FF_IF_Swet'] = Cf_FF_IF_Swet
    results[component] = result
    
    print(f"{component:<12} {Swet:<10.2f} {result['Cf_total']:<8.6f} {FF:<8.3f} {IF:<8.3f} {Cf_FF_IF:<12.6f} {Cf_FF_IF_Swet:<15.6f}")

# Calculate wave drag
CD_wave = wave_drag_coefficient(M, prm.airf_M_crit, prm.MDD)

# Calculate fuselage miscellaneous drag
CD_upsweep = fuselage_upsweep_drag(prm.u_upsweep_deg, prm.A_max_fus, S_w_cur)
CD_base = fuselage_base_drag(M, prm.A_base_fus,  S_w_cur)

# Calculate landing gear drag
CD_landing_gear = landing_gear_drag(gear_deployed, CD_gear)

# Calculate flap drag
CD_flap = flap_drag(flaps_deflected, flap_deflection_deg, c_f, c, S_flap, S_w_cur)

# Summary
print("\n" + "=" * 120)
print("SUMMARY:")
print(f"{'Component':<12} {'Swet [m²]':<10} {'Cf':<8} {'FF':<8} {'IF':<8} {'Cf×FF×IF×Swet':<15}")
print("-" * 70)
for component in components:
    print(f"{component:<12} {results[component]['wetted_area']:<10.2f} {results[component]['Cf_total']:<8.6f} {results[component]['form_factor']:<8.3f} {results[component]['IF']:<8.3f} {results[component]['Cf_FF_IF_Swet']:<15.6f}")


# Calculate drag coefficient if reference area is provided
CD0_friction = total_aircraft_drag_component / S_w_cur

# Add wave drag
print(f"\nWave Drag:")
print(f"Mach number: {M}")
print(f"Wave Drag Coefficient (ΔCD_wave): {CD_wave:.6f}")

# Add fuselage miscellaneous drag
print(f"\nFuselage Miscellaneous Drag:")
print(f"Upsweep angle: {u_upsweep_deg}°")
print(f"Max cross-sectional area (A_max): {A_max_fus} m²")
print(f"Base area (A_base): {A_base_fus} m²")
print(f"Upsweep Drag Coefficient (CD_upsweep): {CD_upsweep:.6f}")
print(f"Base Drag Coefficient (CD_base): {CD_base:.6f}")

# Add landing gear drag
print(f"\nLanding Gear Drag:")
print(f"Landing gear deployed: {'YES' if gear_deployed else 'NO'}")
if gear_deployed:
    print(f"Landing Gear Drag Coefficient (CD_gear): {CD_landing_gear:.6f}")
else:
    print(f"Landing Gear Drag Coefficient: {CD_landing_gear:.6f} (gear retracted)")

# Add flap drag
print(f"\nFlap Drag:")
print(f"Flaps deflected: {'YES' if flaps_deflected else 'NO'}")
if flaps_deflected:
    print(f"Flap deflection: {flap_deflection_deg}°")
    print(f"Flap type: Slotted")
    print(f"Flap chord ratio (c_f/c): {c_f/c:.3f}")
    print(f"Flap area ratio (S_flap/S_w_cur): {S_flap/S_w_cur:.3f}")
    print(f"Flap Drag Coefficient (CD_flap): {CD_flap:.6f}")
else:
    print(f"Flap Drag Coefficient: {CD_flap:.6f} (flaps retracted)")

# Calculate CD0 without excrescence
CD0_without_excrescence = CD0_friction + CD_wave + CD_upsweep + CD_base + CD_landing_gear + CD_flap

# Calculate excrescence and leakage drag
CD_excrescence = excrescence_leakage_drag(CD0_without_excrescence, excrescence_percentage)

# Add excrescence and leakage drag
print(f"\nExcrescence and Leakage Drag:")
print(f"Excrescence percentage: {excrescence_percentage*100:.1f}% of total CD0")
print(f"CD0 without excrescence: {CD0_without_excrescence:.6f}")
print(f"Excrescence & Leakage Drag Coefficient (CD_excrescence): {CD_excrescence:.6f}")

# Total zero-lift drag coefficient
CD0_total = CD0_without_excrescence + CD_excrescence
print(f"  Friction Drag ({CD0_friction:.6f})")
print(f"  + Wave Drag ({CD_wave:.6f})")
print(f"  + Upsweep Drag ({CD_upsweep:.6f})")
print(f"  + Base Drag ({CD_base:.6f})")
print(f"  + Landing Gear Drag ({CD_landing_gear:.6f})")
print(f"  + Flap Drag ({CD_flap:.6f})")
print(f"  + Excrescence & Leakage Drag ({CD_excrescence:.6f})")
print(f"  = Total Zero-Lift Drag Coefficient (CD0): {CD0_total:.6f}")