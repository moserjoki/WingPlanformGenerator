import math
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from sizing import *
from drag import *
from weight import getClassIIWeightEstimation

def C_D0_calculate(flight_condition, gear_deployed, printing):
    # Get atmospheric properties based on flight condition
    atm_properties = get_atmospheric_properties(flight_condition)
    rho = atm_properties['rho']
    altitude = atm_properties['altitude']

    # Flight conditions configuration:
    if flight_condition == 'landing':
        V = prm.V_approach
        M = prm.M_approach
        flow_regime = 'subsonic'
        flaps_deflected = True
        flap_deflection_deg = hld_deflection_land  # Landing flap deflection angle
    elif flight_condition == 'cruise':
        V = prm.V_cruise
        M = prm.M_cruise
        flow_regime = 'transonic'
        flaps_deflected = False
        flap_deflection_deg = 0.0   # No flap deflection in cruise
    elif flight_condition == 'take-off':
        V = math.sqrt(2*m_MTOW*g/(ρ_sea_level*wing.S_w*C_L_max_take_off_cur))
        M = V/343
        flow_regime = 'subsonic'
        flaps_deflected = True
        flap_deflection_deg = hld_deflection_take_off # No flap deflection in cruise

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
        
    # Calculate wave drag
    CD_wave = wave_drag_coefficient(M, prm.airf_M_crit, prm.MDD)

    # Calculate fuselage miscellaneous drag
    CD_upsweep = fuselage_upsweep_drag(prm.u_upsweep_deg, prm.A_max_fus, S_w_cur)
    CD_base = fuselage_base_drag(M, prm.A_base_fus,  S_w_cur)

    # Calculate landing gear drag
    CD_landing_gear = landing_gear_drag(gear_deployed, CD_gear)

    # Calculate flap drag
    CD_flap = flap_drag(flaps_deflected, flap_deflection_deg, c_f, c, S_flap, S_w_cur)

    # Calculate drag coefficient if reference area is provided
    CD0_friction = total_aircraft_drag_component / S_w_cur

    # Calculate CD0 without excrescence
    CD0_without_excrescence = CD0_friction + CD_wave + CD_upsweep + CD_base + CD_landing_gear + CD_flap

    # Calculate excrescence and leakage drag
    CD_excrescence = excrescence_leakage_drag(CD0_without_excrescence, excrescence_percentage)

    # Total zero-lift drag coefficient
    CD0_total = CD0_without_excrescence + CD_excrescence

    if printing:
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

        for component, params in components.items():
            # Calculate skin friction
            print(f"{component:<12} {Swet:<10.2f} {result['Cf_total']:<8.6f} {FF:<8.3f} {IF:<8.3f} {Cf_FF_IF:<12.6f} {Cf_FF_IF_Swet:<15.6f}")

        # Summary
        print("\n" + "=" * 120)
        print("SUMMARY:")
        print(f"{'Component':<12} {'Swet [m²]':<10} {'Cf':<8} {'FF':<8} {'IF':<8} {'Cf×FF×IF×Swet':<15}")
        print("-" * 70)

        for component in components:
            print(f"{component:<12} {results[component]['wetted_area']:<10.2f} {results[component]['Cf_total']:<8.6f} {results[component]['form_factor']:<8.3f} {results[component]['IF']:<8.3f} {results[component]['Cf_FF_IF_Swet']:<15.6f}")

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

        # Add excrescence and leakage drag
        print(f"\nExcrescence and Leakage Drag:")
        print(f"Excrescence percentage: {excrescence_percentage*100:.1f}% of total CD0")
        print(f"CD0 without excrescence: {CD0_without_excrescence:.6f}")
        print(f"Excrescence & Leakage Drag Coefficient (CD_excrescence): {CD_excrescence:.6f}")

        print(f"  Friction Drag ({CD0_friction:.6f})")
        print(f"  + Wave Drag ({CD_wave:.6f})")
        print(f"  + Upsweep Drag ({CD_upsweep:.6f})")
        print(f"  + Base Drag ({CD_base:.6f})")
        print(f"  + Landing Gear Drag ({CD_landing_gear:.6f})")
        print(f"  + Flap Drag ({CD_flap:.6f})")
        print(f"  + Excrescence & Leakage Drag ({CD_excrescence:.6f})")
        print(f"  = Total Zero-Lift Drag Coefficient (CD0): {CD0_total:.6f}")
    return CD0_total

def e_calculate(flap_angle, flaps_deflected, wing_tip_effect, plotting):
    # Drag polar parameters
    CDmin = 0.025
    CLminD = 0.2

    # CL range for the plot
    CL_min = -0.5
    CL_max = 1.5
    # ===================================

    # Calculate aerodynamic factors
    e_factor = oswald_efficiency(wing.AR, wing.quart_sweep, flaps_deflected, flap_angle, wing_tip_effect)
    K_factor = calculate_K(wing.AR, e_factor)

    # Generate CL values for the plot
    CL_values = np.linspace(CL_min, CL_max, 200)

    # Calculate corresponding CD values
    CD_values = drag_polar(CL_values, CDmin, CLminD, K_factor)
    
    if plotting:
        print(f"Calculated parameters:")
        print(f"Oswald efficiency factor (e): {e_factor:.4f}")
        print(f"K factor (1/πAe): {K_factor:.6f}")


        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(CL_values, CD_values, 'b-', linewidth=2, label=f'CD = {CDmin:.3f} + {K_factor:.4f}·(CL - {CLminD:.2f})²')
        plt.axvline(x=CLminD, color='red', linestyle='--', alpha=0.7, label=f'CLminD = {CLminD:.2f}')
        plt.axhline(y=CDmin, color='green', linestyle='--', alpha=0.7, label=f'CDmin = {CDmin:.3f}')

        # Mark the minimum drag point
        min_drag_CL = CLminD
        min_drag_CD = CDmin
        plt.plot(min_drag_CL, min_drag_CD, 'ro', markersize=8, label='Minimum Drag Point')

        # Customize the plot
        plt.xlabel('Lift Coefficient (CL)')
        plt.ylabel('Drag Coefficient (CD)')
        plt.title('Aircraft Drag Polar\n$C_D = C_{D_{min}} + K(C_L - C_{L_{minD}})^2$')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.axis([CL_min, CL_max, 0, max(CD_values) * 1.1])

        # Display the plot
        plt.tight_layout()
        plt.show()
    return e_factor


m_MTOW = 114960 # [kg] Maximum Take Off Weight from Class I Weight Estimation . 
AR = 6.9 # [] aspect ratio

# Create objects used for subsequent calculations
airfoil = Airfoil("airfoils/NASA SC(2)-0414.dat") 
wing = WingSizing(m_MTOW, AR)
cruise_matching_diagram = MatchingDiagram(m_MTOW, AR)

# Dummy values for first matching diagram calculation
C_L_max_take_off_cur = 2.1 # [] Cl during take-off
C_L_max_landing_cur = 2.4 # [] Cl during landing

C_D0_landing_retracted = 0.06 #0.019
C_D0_landing_extended = 0.06 #0.035
C_D0_cruise = 0.06 #0.016
C_D0_take_off_retracted = 0.06# 0.015
C_D0_take_off_extended = 0.06 #0.035

e_landing = 0.8097
e_cruise = 0.626
e_take_off = 0.695

S_w_cur = 0

# Outer loop that takes into account updated m_MTOW and updated C_D0's and e's to rerun aircraft sizing. 
for j in range(6):
    print("\n\n\n")
    print(f"Iteration {j+1}")
    print(f"with m_MTOW: {m_MTOW:03f}")
    wing.update(m_MTOW, AR)
    cruise_matching_diagram.update(m_MTOW, AR)

    # Inner loop that runs the matching diagram to find S_w, then sizes ailerons and HLD's. Runs until C_L_max landing and take off that are required are matched by values that result from calcualtions. 
    printing = False
    for i in range(5):
        if i == 4:
            printing = True

        S_w_cur  = cruise_matching_diagram.compute(C_D0_landing_retracted, C_D0_landing_extended, C_D0_cruise, C_D0_take_off_retracted, C_D0_take_off_extended, 
                                                e_landing, e_cruise, e_take_off, C_L_max_take_off_cur, C_L_max_landing_cur, printing)
        wing.planform_sizing(S_w_cur, printing)
        V_stall = wing.aileron_sizing(C_L_max_landing_cur, printing)

        C_L_max_clean = wing.DATCOM_C_L_max_clean()
        C_L_max_take_off_cur, C_L_max_landing_cur = wing.HLD_sizing(C_L_max_clean)
    
    #cruise_matching_diagram.plot()

    X_cg_aft = 21.98 #RANDOM INITIAL VALUE

    S_v, b_v, c_r_v, c_t_v, MAC_v, Quarter_Chord_Sweep_V, S_h, b_h, c_r_h, c_t_h, MAC_h, Quarter_Chord_Sweep_H = wing.empenage_sizing(X_cg_aft, True)

    C_D0_landing_retracted = C_D0_calculate('landing', False, False)
    C_D0_landing_extended = C_D0_calculate('landing', True, False)
    C_D0_cruise = C_D0_calculate('cruise', False, False)
    C_D0_take_off_retracted = C_D0_calculate('take-off', False, False)
    C_D0_take_off_extended = C_D0_calculate('landing',  True, False)

    e_landing = e_calculate(hld_deflection_land, True, True, False)
    e_cruise = e_calculate(0, False, True, False)
    e_take_off = e_calculate(hld_deflection_take_off, True, True, False)

    print(f"C_D0_landing_retracted: {C_D0_landing_retracted:0.3f} | C_D0_landing_extended {C_D0_landing_extended:0.3f} | C_D0_cruise: {C_D0_cruise:0.3f} | C_D0_take_off_retracted: {C_D0_take_off_retracted:0.3f} | C_D0_take_off_extended: {C_D0_take_off_extended:0.3f}")
    print(f"e_landing {e_landing:0.3f} | e_cruise {e_cruise:0.3f} | e_take_off {e_take_off:0.3f}")

    Ywings = 0.4*l_fus
    Yengine = 0.5*l_fus
    aileronsArea_SI = 4
    subsystem_values = getClassIIWeightEstimation(wing.AR, wing.quart_sweep, wing.taper_ratio, wing.b, wing.S_w, b_h, Ywings, Yengine, S_h, S_v, V_stall, aileronsArea_SI, Quarter_Chord_Sweep_H, Quarter_Chord_Sweep_V)
    
    m_OEW = sum(subsystem_values)*0.453592
    m_payload = 18960 # [kg]
    m_MTOW = (m_OEW + m_payload)*1.45932

    if j == 5:
        wing.fuel_volume(airfoil)
        wing.plot()