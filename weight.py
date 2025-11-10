import math 
import numpy as np
import params as prm
import matplotlib.pyplot as plt








### Constant parameters ###
A_h = prm.empg_AR_h       # ???? Horizontal tail Aspect ratio, [-]
A_v = prm.empg_AR_v     # Vertical Tail Aspect ratio, [-]

D = 13.559        # fuselage structural depth, [ft]
L =  145.906329       # fuselage structural length (excludes radome, tail cap), [ft]

K_door = 1.06       # For one side cargo door
K_Lg = 1.0      #1.12 if fuselage-mounted main landing gear; =1.0 otherwise
K_ng = 1.017        # for pylon-mounted nacelle
K_np, K_mp = 1.0, 1.0      #1.126 for kneeling gear; =1.0 otherwise
K_p = 1.0       # 1.4 for engine with propeller
K_r = 1.0       # 1.133 for reciprocating engine 
K_tp = 1.0      # 0.793 if turboprop
K_tr = 1.18     # 1.18 for engine with thrust reverser, 1.0 otherwise
K_uht = 1.0 # 1.143 for unit (all-moving) horizontal tail; 1.0 otherwise
H_t_over_H_v = 0        #0.0 for conventional tail, [-]
N_c = 11      # number of crew, [-]
N_en = 2         # number of engines, [-]
N_f = 7     # number of functions performed by controls (typically 4-7)
N_l =       # ultimate landing load factor; =N_gear*1.5 [-]
N_Lt = 11.7833333     # nacelle length, [ft]
#The documentation provides flange to flange length as 141.4 inches = 11.7833333 ft
#https://prd-sc102-cdn.rtx.com/-/media/pw/products/commercial-jet-engines/pw2000/files/ce_pw2000_fact.pdf?rev=-1&hash=3E2400E9D1BAA5B7B22E700214321D31 

N_mss = 2    # number of main gear shock struts, [-]
N_m = 0     # number of mechanical functions (typically 0-2)
N_mw = 4     # number of main wheels, [-]
N_nw = 2     # number of nose wheels, [-]
N_p = 188 + N_c      # number of personnel onboard (crew and passagers)

N_t =  3     # number of fuel tanks, [-]
#one in each wing and one inside the fuselage

N_w = 6.54166667      # nacelle width, [ft]
#The documentation provides fan tip diameter as 78.5 inches = 6.541667 feet
#https://prd-sc102-cdn.rtx.com/-/media/pw/products/commercial-jet-engines/pw2000/files/ce_pw2000_fact.pdf?rev=-1&hash=3E2400E9D1BAA5B7B22E700214321D31

S_cs =        # control surface area (all), [ft^2]
S_n = 242.162246       # nacelle wetted area, [ft^2]
#Assumed surface area of a cylinder with length N_Lt and diameter N_w therefore wetted surface area is given by 2 pi r h

R_kva = 55  # system electrocal rating, typical values for cargo aircrafts
t_c_root = 0.122     # based on chosen airfoil
V_i = 11346.19      # integral tanks volume, [gal]
V_t = 11346.19      # total fuel volume, [gal]
V_p = 0     # self-sealing "protected" tanks volume, [gal], apparently only military aircraft
W_APUUninstalled = 280 # [lb]
# Honeywell HGT1700, APU used in Airbus A350
W_c = 8289.38106     # Maximum cargo weight, [lb]
#The average mass per passenger including luggage was given as 98.8kg, I assumed 20kg of it to be the luggage(cargo) mass per passenger
#Therefore total maximum cargo weight is 188 * 20kg = 8289lbs
W_en = 7299.946425     # engine weight, each, [lb]
W_fw = 75750.833     # weight of fuel in wing, [lb]
W_uav = 1200    # uninstalled avionics weight, [lb]



def getClassIIWeightEstimation(
        aspectRatio,
        sweepWings,
        taperWings,
        wingSpan_SI,
        wingArea_SI,

        horizontalTailSpan_SI,

        YhorizontalTail,
        Ywings,
        Yengine
        ):
### Input parameters ###
    A=  aspectRatio  # aspect ratio [-]
    B_w = wingSpan_SI*3.2808399     # wing span, [ft]
    B_h = horizontalTailSpan_SI*3.2808399    # horizontal tail span, [ft]
    
    L_f = prm.l_fus*3.2808399      # total fuselage length, [ft] ????
    L_t = (YhorizontalTail-Ywings)*3.2808399      # tail length, wing quater-MAC to tail-quater-MAC, [ft]
    L_m =       # length of main landing gear, [in]
    L_n =       # length of nose landing gear, [in]
    N_z=        # ultimate load factor; 1.5* limit load factor, [-]
    S_csw =     # control surface area (wing-mounted), [ft^2]

    S_f = (prm.l_fus)*(prm.d_fuselage)*np.pi*10.7639      # Fuselage wetted area, [ft^2]
    S_ht =      # Horizontal tail area, [ft^2]
    S_vt =      # Vertical tail area, [ft^2]
    S_w= wingArea_SI*10.7639104      #trapezoidal wing area, [ft^2]

    V_stall =       # Stall speed, [ft/s]?????

    W_l = 0.84*253443      # Landing design gross weight, [lb]

    Lambda = sweepWings      # wing sweep at 1/4 (25%) MAC [deg]
    lambda_ = taperWings     # Taper ratio [-]
    Lambda_ht =  max(34.2, sweepWings)      # horizontal tail sweep at 1/4 (25%) MAC [deg]
    Lambda_vt = max(37.8)# horizontal tail sweep at 1/4 (25%) MAC [deg]



    ### Parameters based on input parameters ###
    F_w = (prm.l_fus-YhorizontalTail)/(prm.L3_fus)*prm.d_fuselage*3.2808399 # Fuselage width at horizontal tail intersection, [ft]
    K_y = 0.3*L_t      #aircraft pitching radius of gyration, [ft] (approx 0.3 L_t)
    K_z = 1*L_t        # Aircraft Yawing radius of gyration, [ft] (approx L_t)
    I_y = W_dg*K_z**2 # yawing moment of inertia [lb ft^2]
    #calculation based on the approximate radius of gyration
    L_ec = Yengine * 3.2808399      # length from engine front to cockpit (total if multiengine), [ft]
    S_e =       #elevator area, [ft^2]
    W_ec = 2.331*W_en**0.901*K_p*K_tr     # weight of engine and contents, (per nacelle), [lb]
    K_ws = 0.75*((1+2*lambda_)/(1+lambda_))*(B_w * np.tan(np.radians(Lambda))/L)
    N_gen = N_en    # number of generators (typically =N_en)
    L_a =   2*(L_ec+0.35*0.5*B_w)    # electrical routing distance, generators to avionics to cockpit, [ft]
    #approximated as the distance from engines to centerline and to cockpit 
    V_pr =  0.95*(prm.l_fus)*(prm.d_fuselage/2)**2*np.pi*35.3146667     # volume of pressurized section [ft^3]
    #approximated as the volume of a cylinder with length equal to the fuselage length and diameter equal to the fuselage diameter, with 0.95

    W_dg= 253443   # Design gross weight, [lb]
    # Result of Class I weight estimation



    ### Functions ###
    W_wing = 0.0051*(W_dg*N_z)**0.557 * S_w**0.649 * A**0.5 * t_c_root**(-0.4) * (1+lambda_)**0.1 * np.cos(np.radians(Lambda))**(-1)*S_csw**0.1

    W_horizontalTail = 0.0379 * K_uht * (1+F_w/B_h)**(-0.25) * W_dg**0.639 * N_z**0.1 * S_ht**0.75 *L_t**(-1) * K_y**0.704 * (np.cos(np.radians(Lambda_ht)))**(-1) * A_h**0.166 * (1+S_e/S_ht)**0.1

    W_verticalTail = 0.0026 * (1+H_t_over_H_v)**0.225 * W_dg**0.556 * N_z**536 * L_t**(-0.5) * S_vt**0.5 * K_z**0.875 * np.cos(np.radians(Lambda_vt))**(-1) * A_v**0.35 * (t_c_root)**(-0.5)

    W_fuselage = 0.3280 * K_door * K_Lg *(W_dg*N_z)**0.5 * L**0.25 * S_f**0.302 * (1+K_ws)**0.04 * (L/D)**0.10

    W_mainLandingGear = 0.0106* K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**(-0.5) * V_stall**0.1

    W_noseLandingGear = 0.032 * K_np * W_l**0.646 * N_l**0.2 * L_n**0.5 * N_nw **0.45

    W_nacelleGroup = 0.6724 * K_ng * N_Lt**0.10 * N_w**0.294 * N_z**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224

    W_engineControls = 5.0*N_en + 0.80*L_ec

    W_starterPneumatic = 49.19*(N_en*W_en/1000)**0.541

    W_fuelSystem = 2.405 * V_t**0.606 * (1+V_i/V_t)*(-1) * (1 + V_p/V_t)*N_t**0.5

    W_flightControls = 145.9 * N_f**0.554 * (1 + N_m/N_f)**(-1) * S_cs**0.20 * (I_y * 10**(-6))*0.07

    W_APUInstalled = 2.2 * W_APUUninstalled 

    W_instruments = 4.509 * K_r * K_tp * N_c**0.541 * N_en * (L_f + B_w)**0.5

    W_hydraulics = 0.2673 * 0.2673 * N_f * (L_f + B_w)**0.937

    W_electrical = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10

    W_avionics = 1.73 * W_uav**0.983

    W_furnishings = 0.0577*N_c**0.1*W_c**0.393*S_f**0.75

    W_airConditioning = 62.36*N_p**0.25*(V_pr/1000)**0.604*W_uav**0.10

    W_antiIce = 0.002*W_dg

    W_handling_gear = 0.0003 * W_dg 
    
    # Ordered subsystem names and current-iteration values (list approach)

    subsystem_values = [
        W_wing,
        W_horizontalTail,
        W_verticalTail,
        W_fuselage,
        W_mainLandingGear,
        W_noseLandingGear,
        W_nacelleGroup,
        W_engineControls,
        W_starterPneumatic,
        W_fuelSystem,
        W_flightControls,
        W_APUInstalled,
        W_instruments,
        W_hydraulics,
        W_electrical,
        W_avionics,
        W_furnishings,
        W_airConditioning,
        W_antiIce,
        W_handling_gear,
    ]
    return subsystem_values


def plotWeightBreakdown(subsystem_names, subsystem_values):
    subsystem_names = [
                "Wing",
                "Horizontal Tail",
                "Vertical Tail",
                "Fuselage",
                "Main Landing Gear",
                "Nose Landing Gear",
                "Nacelle Group",
                "Engine Controls",
                "Starter Pneumatic",
                "Fuel System",
                "Flight Controls",
                "APU Installed",
                "Instruments",
                "Hydraulics",
                "Electrical",
                "Avionics",
                "Furnishings",
                "Air Conditioning",
                "Anti-Ice",
                "Handling Gear",
    ]
    #Create a bar chart to visualize the weight breakdown
    plt.figure(figsize=(12, 6))
    plt.bar(subsystem_names, subsystem_values)
    plt.xlabel('Subsystems')
    plt.ylabel('Weight (lb)')
    plt.title('Aircraft Weight Breakdown by Subsystem')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    