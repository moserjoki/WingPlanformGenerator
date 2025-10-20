#ZERO-LIFT DRAG COEFFICIENT

from math import *
import params as prm

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
def Swet_wing(b, c_tip, c_fuselage_intersect, d_fus):
    S_expW = ((b-d_fus)/2)*(c_tip+c_fuselage_intersect)
    S_wetW = 1.07 * 2 * S_expW
    return S_wetW 

def Swet_HT(b_ht, c_tip_ht, c_root_ht):
    S_expHT = (b_ht/2)*((c_tip_ht+c_root_ht)/2)
    S_wetHT = 1.05 * 2 * S_expHT
    return S_wetHT

def Swet_VT(b_vt, c_tip_vt, c_root_vt):
    S_expVT = (b_vt)*((c_tip_vt+c_root_vt)/2)
    S_wetVT = 1.05 * 2 * S_expVT
    return S_wetVT

def Swet_fus():
    term1 = (1 / (3 * prm.L1_fus**2)) * ((4 * prm.L1_fus**2 + prm.d_fuselage**2/4)**1.5-prm.d_fuselage**3/8)
    term2 = -prm.d_fuselage + 4 * prm.L2_fus + 2 * sqrt(prm.L3_fus**2 + prm.d_fuselage**2/4)
    Swet = (pi*prm.d_fuselage/4) * (term1 + term2)
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


