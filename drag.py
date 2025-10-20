#ZERO-LIFT DRAG COEFFICIENT

from math import *

# === Atmospheric Properties ===
def get_atmospheric_properties(flight_condition):
    """
    Get air density and dynamic viscosity for different flight conditions.
    
    Parameters:
    flight_condition (str): 'takeoff' or 'cruise'
    
    Returns:
    dict: Atmospheric properties
    """
    if flight_condition == 'takeoff':
        # Sea level conditions (ISA)
        return {
            'rho': 1.225,      # Air density [kg/m³]
            'mu': 1.789e-5,    # Dynamic viscosity [kg/(m·s)]
            'altitude': 'Sea Level'
        }
    elif flight_condition == 'cruise':
        # Cruise at 31,000 ft (approximately 9,450 meters)
        return {
            'rho': 0.459,      # Air density [kg/m³] at 31,000 ft
            'mu': 1.427e-5,    # Dynamic viscosity [kg/(m·s)] at 31,000 ft
            'altitude': '31,000 ft (9,450 m)'
        }
    else:
        # Default to takeoff conditions
        return {
            'rho': 1.225,
            'mu': 1.789e-5,
            'altitude': 'Sea Level'
        }

# === Wetted Areas ===
def Swet_wing(S_expW):
    return 1.07 * 2 * S_expW

def Swet_HT(S_expHT):
    return 1.05 * 2 * S_expHT

def Swet_VT(S_expVT):
    return 1.05 * 2 * S_expVT

def Swet_fus(d_fus, L1, L2, L3):
    term1 = (1 / (3 * L1**2)) * ((4 * L1**2 + d_fus**2 / 4)**1.5 - d_fus**3 / 8)
    term2 = -d_fus + 4 * L2 + 2 * sqrt(L3**2 + d_fus**2 / 4)
    Swet = (pi * d_fus / 4) * (term1 + term2)
    return Swet

# === Skin friction coefficients ===

def calculate_component_skin_friction(rho, V, l, k, M, flow_regime, laminar_fraction, flight_condition='takeoff'):
    """
    Calculate skin friction coefficient for a single component.
    """
    # Get atmospheric properties based on flight condition
    atm_properties = get_atmospheric_properties(flight_condition)
    mu = atm_properties['mu']
    
    # Calculate actual Reynolds number
    Re_actual = (rho * V * l) / mu
    
    # Calculate cutoff Reynolds number based on flow regime
    if flow_regime == 'subsonic':
        Re_cutoff = 38.21 * (l / k) ** 1.053
    else:  # transonic
        Re_cutoff = 44.62 * (l / k) ** 1.053 * M ** 1.16
    
    # Take minimum value
    Re_effective = min(Re_actual, Re_cutoff)
    
    # Calculate laminar and turbulent skin friction coefficients
    Cf_laminar = 1.328 / Re_effective
    
    log_re = log10(Re_effective)
    Cf_turbulent = 0.455 / (log_re ** 2.58)
    
    # Apply Mach correction for turbulent flow if M > 0.3
    if M > 0.3:
        Cf_turbulent /= (1 + 0.144 * M ** 2) ** 0.65
    
    # Calculate weighted average
    Cf_total = (laminar_fraction * Cf_laminar) + ((1 - laminar_fraction) * Cf_turbulent)
    
    return {
        'Cf_total': Cf_total,
        'Cf_laminar': Cf_laminar,
        'Cf_turbulent': Cf_turbulent,
        'Re_actual': Re_actual,
        'Re_cutoff': Re_cutoff,
        'Re_effective': Re_effective
    }

# === Form Factor Functions ===

def form_factor_lifting_surfaces(t_c, x_c_max, M, Lambda_max_deg):
    """
    Form factor for lifting surfaces (wing, tail, strut, pylon)
    
    Parameters:
    t_c (float): Thickness to chord ratio
    x_c_max (float): Location of maximum thickness (x/c)
    M (float): Mach number
    Lambda_max_deg (float): Sweep angle at maximum thickness line [degrees]
    """
    # Convert sweep angle to radians for cosine
    Lambda_max_rad = radians(Lambda_max_deg)
    
    # Basic form factor
    FF_basic = 1 + (0.6 / x_c_max) * t_c + 100 * (t_c ** 4)
    
    # Compressibility correction (only for M > 0.2)
    if M > 0.2:
        FF_compressibility = 1.34 * (M ** 0.18) * (cos(Lambda_max_rad) ** 0.28)
        FF_total = FF_basic * FF_compressibility
    else:
        FF_total = FF_basic
    
    return FF_total

def form_factor_fuselage(l_fuselage, A_max, Re_effective):
    """
    Form factor for fuselage and smooth canopy
    
    Parameters:
    l_fuselage (float): Length of fuselage [m]
    A_max (float): Maximum cross-sectional (frontal) area [m²]
    Re_effective (float): Effective Reynolds number
    """
    # Only apply for supercritical Reynolds > 10^5
    if Re_effective > 1e5:
        # Calculate fineness ratio
        d_equivalent = 2 * sqrt(A_max / pi)  # Equivalent diameter
        f = l_fuselage / d_equivalent
        
        FF = 1 + (60 / (f ** 3)) + (f / 400)
    else:
        FF = 1.0  # Default for subcritical Reynolds
    
    return FF

def form_factor_nacelle(l_nacelle, A_max):
    """
    Form factor for nacelle and smooth external store
    
    Parameters:
    l_nacelle (float): Length of nacelle [m]
    A_max (float): Maximum cross-sectional area [m²]
    """
    # Calculate fineness ratio
    d_equivalent = 2 * sqrt(A_max / pi)  # Equivalent diameter
    f = l_nacelle / d_equivalent
    
    FF = 1 + (0.35 / f)
    
    return FF

# === Wave Drag Function ===
def wave_drag_coefficient(M, M_cr, M_DD):
    """
    Calculate wave drag coefficient based on Mach number.
    
    Parameters:
    M (float): Current Mach number
    M_cr (float): Critical Mach number
    M_DD (float): Drag divergence Mach number
    
    Returns:
    float: Wave drag coefficient increment
    """
    if M < M_cr:
        # No wave drag below critical Mach
        return 0.0
    elif M_cr <= M <= M_DD:
        # Between critical and drag divergence Mach
        return 0.002 * (1 + 2.5 * (M_DD - M) / 0.05) ** -1
    else:
        # Above drag divergence Mach
        return 0.002 * (1 + (M - M_DD) / 0.05) ** 2.5

# === Fuselage Miscellaneous Drag Functions ===
def fuselage_upsweep_drag(u_deg, A_max, S_ref):
    """
    Calculate fuselage upsweep drag coefficient.
    
    Parameters:
    u_deg (float): Upsweep angle [degrees]
    A_max (float): Maximum fuselage cross-sectional area [m²]
    S_ref (float): Reference wing area [m²]
    
    Returns:
    float: Upsweep drag coefficient
    """
    # Convert upsweep angle from degrees to radians
    u_rad = radians(u_deg)
    
    # Calculate D/q
    D_q = 3.83 * (u_rad ** 2.5) * A_max
    
    # Convert to coefficient by dividing by S_ref
    CD_upsweep = D_q / S_ref
    
    return CD_upsweep

def fuselage_base_drag(M, A_base, S_ref):
    """
    Calculate fuselage base drag coefficient (subsonic only).
    
    Parameters:
    M (float): Mach number
    A_base (float): Fuselage base area [m²]
    S_ref (float): Reference wing area [m²]
    
    Returns:
    float: Base drag coefficient
    """
    # Only for subsonic Mach numbers
    if M >= 1.0:
        return 0.0  # Formula only valid for subsonic
    
    # Calculate D/q
    D_q = 0.139 + 0.419 * ((M - 0.161) ** 2)
    D_q *= A_base
    
    # Convert to coefficient by dividing by S_ref
    CD_base = D_q / S_ref
    
    return CD_base

# === Landing Gear Drag Function ===
def landing_gear_drag(gear_deployed, CD_gear_constant=0.0):
    """
    Calculate landing gear drag contribution.
    
    Parameters:
    gear_deployed (bool): Whether landing gear is deployed
    CD_gear_constant (float): Constant drag coefficient for landing gear
    
    Returns:
    float: Landing gear drag coefficient
    """
    if gear_deployed:
        return CD_gear_constant
    else:
        return 0.0

# === Flap Drag Function ===
def flap_drag(flaps_deflected, flap_deflection_deg, c_f, c, S_flap, S_ref):
    """
    Calculate flap drag coefficient (slotted flaps only).
    
    Parameters:
    flaps_deflected (bool): Whether flaps are deflected
    flap_deflection_deg (float): Flap deflection angle [degrees]
    c_f (float): Chord length of flap [m]
    c (float): Wing chord length [m]
    S_flap (float): Wing area affected by flap [m²]
    S_ref (float): Reference wing area [m²]
    
    Returns:
    float: Flap drag coefficient
    """
    if not flaps_deflected or flap_deflection_deg <= 10:
        return 0.0  # No flap drag if not deflected or deflection ≤ 10°
    
    # Slotted flap constant
    F_flap = 0.0074
    
    # Calculate flap drag coefficient
    CD_flap = F_flap * (c_f / c) * (S_flap / S_ref) * (flap_deflection_deg - 10)
    
    return CD_flap

# === Excrescence and Leakage Drag Function ===
def excrescence_leakage_drag(CD0_without_excrescence, excrescence_percentage=0.035):
    """
    Calculate excrescence and leakage drag as percentage of total CD0.
    
    Parameters:
    CD0_without_excrescence (float): CD0 without excrescence drag
    excrescence_percentage (float): Percentage as decimal (0.02-0.05 for 2-5%)
    
    Returns:
    float: Excrescence and leakage drag coefficient
    """
    return CD0_without_excrescence * excrescence_percentage

# === INPUT PARAMETERS FOR EACH COMPONENT ===

# Flight conditions (same for all components)
flight_condition = 'takeoff'  # 'takeoff' or 'cruise'

# Get atmospheric properties based on flight condition
atm_properties = get_atmospheric_properties(flight_condition)
rho = atm_properties['rho']
altitude = atm_properties['altitude']

# Set appropriate velocities and Mach numbers based on flight condition
if flight_condition == 'takeoff':
    V = 80.0          # Takeoff velocity [m/s]
    M = 0.23          # Takeoff Mach number
    flow_regime = 'subsonic'
    gear_deployed = True
    flaps_deflected = True
    flap_deflection_deg = 20.0  # Takeoff flap setting
else:  # cruise
    V = 230.0         # Cruise velocity [m/s]
    M = 0.78          # Cruise Mach number
    flow_regime = 'transonic'
    gear_deployed = False
    flaps_deflected = False
    flap_deflection_deg = 0.0   # No flap deflection in cruise

k = 1e-5         # Skin roughness [m] (same for all components)

# Wave drag parameters
M_cr = 0.75      # Critical Mach number
M_DD = 0.82      # Drag divergence Mach number

# Fuselage miscellaneous drag parameters
u_upsweep_deg = 5.0     # Fuselage upsweep angle [degrees]
A_max_fus = 2.5         # Maximum fuselage cross-sectional area [m²]
A_base_fus = 0.8        # Fuselage base area [m²]

# Landing gear parameters
CD_gear = 0.015         # Constant drag coefficient for landing gear

# Flap parameters (slotted flaps)
c_f = 0.8               # Chord length of flap [m]
c = 2.0                 # Wing chord length [m]
S_flap = 10.0           # Wing area affected by flap [m²]

# Excrescence and leakage parameters
excrescence_percentage = 0.035  # 3.5% of total CD0 (2-5% range)

# Geometry parameters for wetted areas
S_expW = 15.0    # Exposed wing area [m²]
S_expHT = 3.0    # Exposed horizontal tail area [m²]
S_expVT = 2.0    # Exposed vertical tail area [m²]
d_fus = 1.8      # Fuselage diameter [m]
L1 = 5.0         # Forward fuselage length [m]
L2 = 8.0         # Midsection length [m]
L3 = 2.0         # Aft fuselage length [m]

# Reference area
S_ref = 15.0     # Reference wing area [m²]

# Calculate wetted areas
S_wet_wing = Swet_wing(S_expW)
S_wet_HT = Swet_HT(S_expHT)
S_wet_VT = Swet_VT(S_expVT)
S_wet_fus = Swet_fus(d_fus, L1, L2, L3)

# Component-specific parameters with interference factors and wetted areas
components = {
    'wing': {
        'l': 2.0,           # Characteristic length [m] (mean aerodynamic chord)
        'laminar_fraction': 0.1,   # 10% laminar flow on wing
        'type': 'lifting_surface',
        't_c': 0.12,        # Thickness to chord ratio
        'x_c_max': 0.3,     # Location of max thickness
        'Lambda_max_deg': 25.0,  # Sweep at max thickness [deg]
        'IF': 1.05,  # Wing-fuselage interference
        'wetted_area': S_wet_wing
    },
    'fuselage': {
        'l': 15.0,          # Fuselage length [m]
        'laminar_fraction': 0.02,  # 2% laminar flow on fuselage
        'type': 'fuselage',
        'A_max': A_max_fus, # Maximum cross-sectional area [m²]
        'IF': 1.10,  # Fuselage base drag, protuberances
        'wetted_area': S_wet_fus
    },
    'HT': {
        'l': 1.2,           # Horizontal tail chord [m]
        'laminar_fraction': 0.05,  # 5% laminar flow on HT
        'type': 'lifting_surface', 
        't_c': 0.10,        # Thickness to chord ratio
        'x_c_max': 0.3,     # Location of max thickness
        'Lambda_max_deg': 15.0,  # Sweep at max thickness [deg]
        'IF': 1.10,  # Tail-fuselage interference
        'wetted_area': S_wet_HT
    },
    'VT': {
        'l': 1.5,           # Vertical tail chord [m]
        'laminar_fraction': 0.05,  # 5% laminar flow on VT
        'type': 'lifting_surface',
        't_c': 0.10,        # Thickness to chord ratio
        'x_c_max': 0.3,     # Location of max thickness
        'Lambda_max_deg': 35.0,  # Sweep at max thickness [deg]
        'IF': 1.10,  # Tail-fuselage interference
        'wetted_area': S_wet_VT
    }
}

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
CD_wave = wave_drag_coefficient(M, M_cr, M_DD)

# Calculate fuselage miscellaneous drag
CD_upsweep = fuselage_upsweep_drag(u_upsweep_deg, A_max_fus, S_ref)
CD_base = fuselage_base_drag(M, A_base_fus, S_ref)

# Calculate landing gear drag
CD_landing_gear = landing_gear_drag(gear_deployed, CD_gear)

# Calculate flap drag
CD_flap = flap_drag(flaps_deflected, flap_deflection_deg, c_f, c, S_flap, S_ref)

# Summary
print("\n" + "=" * 120)
print("SUMMARY:")
print(f"{'Component':<12} {'Swet [m²]':<10} {'Cf':<8} {'FF':<8} {'IF':<8} {'Cf×FF×IF×Swet':<15}")
print("-" * 70)
for component in components:
    print(f"{component:<12} {results[component]['wetted_area']:<10.2f} {results[component]['Cf_total']:<8.6f} {results[component]['form_factor']:<8.3f} {results[component]['IF']:<8.3f} {results[component]['Cf_FF_IF_Swet']:<15.6f}")

print(f"\nTotal Aircraft Drag Component (sum of Cf × FF × IF × Swet): {total_aircraft_drag_component:.6f} m²")

# Calculate equivalent flat plate drag area
print(f"Equivalent Flat Plate Drag Area: {total_aircraft_drag_component:.6f} m²")

# Calculate drag coefficient if reference area is provided
CD0_friction = total_aircraft_drag_component / S_ref
print(f"Friction Drag Coefficient (based on S_ref = {S_ref} m²): {CD0_friction:.6f}")

# Add wave drag
print(f"\nWave Drag:")
print(f"Mach number: {M}")
print(f"Critical Mach (M_cr): {M_cr}")
print(f"Drag Divergence Mach (M_DD): {M_DD}")
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
    print(f"Flap area ratio (S_flap/S_ref): {S_flap/S_ref:.3f}")
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
print(f"\nTotal Zero-Lift Drag Coefficient (CD0): {CD0_total:.6f}")
print(f"  = Friction Drag ({CD0_friction:.6f})")
print(f"  + Wave Drag ({CD_wave:.6f})")
print(f"  + Upsweep Drag ({CD_upsweep:.6f})")
print(f"  + Base Drag ({CD_base:.6f})")
print(f"  + Landing Gear Drag ({CD_landing_gear:.6f})")
print(f"  + Flap Drag ({CD_flap:.6f})")
print(f"  + Excrescence & Leakage Drag ({CD_excrescence:.6f})")