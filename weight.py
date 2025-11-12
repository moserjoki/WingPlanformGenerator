import math 
import numpy as np
import params as prm
import matplotlib.pyplot as plt

def getClassIIWeightEstimation(
        aspectRatio,
        sweepWings,
        taperWings,
        wingSpan_SI,
        wingArea_SI,

        horizontalTailSpan_SI,

        Ywings,
        Yengine,
        HorizontalTail_area,
        VerticalTail_area,
        V_Stall,

        aileronsArea_SI, 
        Quarter_Chord_Sweep_H, 
        Quarter_Chord_Sweep_V,

        MTOW_initial_SI
        ):
    ### Constant parameters ###
    A_h = prm.empg_AR_h       # ???? Horizontal tail Aspect ratio, [-]
    A_v = prm.empg_AR_v     # Vertical Tail Aspect ratio, [-]

    D = 13.559        # fuselage structural depth, [ft]
    L = 145.906329    # fuselage structural length (excludes radome, tail cap), [ft]

    K_door = 1.06   # For one side cargo door
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
    N_l =  1.8     # ultimate landing load factor; =N_gear*1.5 [-]
    # estimated based on Comet and EASA
    N_Lt = 11.7833333     # nacelle length, [ft]
    #The documentation provides flange to flange length as 141.4 inches = 11.7833333 ft
    #https://prd-sc102-cdn.rtx.com/-/media/pw/products/commercial-jet-engines/pw2000/files/ce_pw2000_fact.pdf?rev=-1&hash=3E2400E9D1BAA5B7B22E700214321D31 

    N_mss = 2    # number of main gear shock struts, [-]
    N_m = 0     # number of mechanical functions (typically 0-2)
    N_mw = 4     # number of main wheels, [-]
    N_nw = 2     # number of nose wheels, [-]
    N_p = 188 + N_c      # number of personnel onboard (crew and passagers)

    N_t =  4     # number of fuel tanks, [-]
    #one in each wing and one inside the fuselage

    N_w = 6.54166667      # nacelle width, [ft]
    #The documentation provides fan tip diameter as 78.5 inches = 6.541667 feet
    #https://prd-sc102-cdn.rtx.com/-/media/pw/products/commercial-jet-engines/pw2000/files/ce_pw2000_fact.pdf?rev=-1&hash=3E2400E9D1BAA5B7B22E700214321D31


    S_n = 242.162246       # nacelle wetted area, [ft^2]
    #Assumed surface area of a cylinder with length N_Lt and diameter N_w therefore wetted surface area is given by 2 pi r h

    R_kva = 55  # system electrocal rating, typical values for cargo aircrafts
    t_c_root = 0.122     # based on chosen airfoil

    V_p = 0     # self-sealing "protected" tanks volume, [gal], apparently only military aircraft
    W_APUUninstalled = 280 # [lb]
    # Honeywell HGT1700, APU used in Airbus A350
    W_c = 8289.38106     # Maximum cargo weight, [lb]
    #The average mass per passenger including luggage was given as 98.8kg, I assumed 20kg of it to be the luggage(cargo) mass per passenger
    #Therefore total maximum cargo weight is 188 * 20kg = 8289lbs
    W_en = 7300     # engine weight, each, [lb]

    W_uav = 1200    # uninstalled avionics weight, [lb]


### Input parameters ###
    fuel_volume_initial_SI = (MTOW_initial_SI/1.459328551*0.459328551)/800
    A=  aspectRatio  # aspect ratio [-]
    B_w = wingSpan_SI*3.2808399     # wing span, [ft]
    B_h = horizontalTailSpan_SI*3.2808399    # horizontal tail span, [ft]
    
    L_f = prm.l_fus*3.2808399      # total fuselage length, [ft] ????
    YhorizontalTail = 0.9*prm.l_fus #  Distance from the nose to the quarter chord point of the MAC of the horizontal tail [m]
    L_t = (YhorizontalTail-Ywings)*3.2808399      # tail length, wing quater-MAC to tail-quater-MAC, [ft]
    L_m =  138     # length of main landing gear, [in]
    L_n =  42    # length of nose landing gear, [in]
    ### achtung! this is reasonable estimage based on Comet, need to be updated with better data ###

    N_z =  3.75     # ultimate load factor; 1.5* limit load factor, [-]
    
    

    S_f = 4848.17289     # Fuselage wetted area, [ft^2]
    S_ht =  HorizontalTail_area*10.7639104   # Horizontal tail area, [ft^2]
    S_vt =  VerticalTail_area*10.7639104  # Vertical tail area, [ft^2]
    S_w= wingArea_SI*10.7639104      #trapezoidal wing area, [ft^2]
    print("Wing area ft2:", S_w, "Horizontal tail area ft2:", S_ht, "Vertical tail area ft2:", S_vt)

    V_stall =    V_Stall *3.28084  # Stall speed, [ft/s]

    V_i = fuel_volume_initial_SI*264.172052    # integral tanks volume, [gal]
    V_t = V_i   # total fuel volume, [gal]

    W_l = 0.84*MTOW_initial_SI*2.20462262   # Landing design gross weight, [lb]

    Lambda = sweepWings      # wing sweep at 1/4 (25%) MAC [deg]
    lambda_ = taperWings     # Taper ratio [-]
    Lambda_ht = Quarter_Chord_Sweep_H # horizontal tail sweep at 1/4 (25%) MAC [deg]
    Lambda_vt = Quarter_Chord_Sweep_V # vertical tail sweep at 1/4 (25%) MAC [deg]



    ### Parameters based on input parameters ###
    F_w = (prm.l_fus-YhorizontalTail)/(prm.L3_fus)*prm.d_fuselage*3.2808399 # Fuselage width at horizontal tail intersection, [ft]
    K_y = 0.3*L_t      #aircraft pitching radius of gyration, [ft] (approx 0.3 L_t)
    K_z = 1*L_t        # Aircraft Yawing radius of gyration, [ft] (approx L_t)
    W_dg= MTOW_initial_SI*2.20462262   # Design gross weight, [lb]
    I_y = 1800000 # yawing moment of inertia [lb ft^2]
    #based on reference aircraft data given by Comet
    L_ec = Yengine * 1.5 * 2      # length from engine front to cockpit (total if multiengine), [ft]
    S_e = 0.25*S_ht      #elevator area, [ft^2]
    #approximated as 25% of horizontal tail area
    W_ec = 2.331*W_en**0.901*K_p*K_tr     # weight of engine and contents, (per nacelle), [lb]
    K_ws = 0.75*((1+2*lambda_)/(1+lambda_))*(B_w * np.tan(np.radians(Lambda))/L)
    N_gen = N_en    # number of generators (typically =N_en)
    L_a =   2*(L_ec+0.35*0.5*B_w)    # electrical routing distance, generators to avionics to cockpit, [ft]
    #approximated as the distance from engines to centerline and to cockpit 
    V_pr =  0.95*(prm.l_fus)*(prm.d_fuselage/2)**2*np.pi*35.3146667     # volume of pressurized section [ft^3]
    #approximated as the volume of a cylinder with length equal to the fuselage length and diameter equal to the fuselage diameter, with 0.95

    W_fw = fuel_volume_initial_SI*800*2.20462262   # weight of fuel in wing, [lb]
    #Only if all fuel can be stored in the wings, should be checked

    # Result of Class I weight estimation

    S_csw =  aileronsArea_SI*10.7639104   # control surface area (wing-mounted), [ft^2]
    S_cs = S_e + S_csw +  HorizontalTail_area*0.25*10.7639104  # control surface area (all), [ft^2]


    ### Functions ###
    W_wing = 0.0051*(W_dg*N_z)**0.557 * S_w**0.649 * A**0.5 * t_c_root**(-0.4) * (1+lambda_)**0.1 * np.cos(np.deg2rad(Lambda))**(-1)*S_csw**0.1

    W_horizontalTail = 0.0379 * K_uht * (1+F_w/B_h)**(-0.25) * W_dg**0.639 * N_z**0.1 * S_ht**0.75 *L_t**(-1) * K_y**0.704 * (np.cos(np.deg2rad(Lambda_ht)))**(-1) * A_h**0.166 * (1+S_e/S_ht)**0.1

    W_verticalTail = 0.0026 * (1+H_t_over_H_v)**0.225 * W_dg**(0.556) * N_z**0.536 * L_t**(-0.5) * S_vt**0.5 * K_z**0.875 * np.cos(np.radians(Lambda_vt))**(-1) * A_v**0.35 * (t_c_root)**(-0.5)

    W_fuselage = 0.3280 * K_door * K_Lg *(W_dg*N_z)**0.5 * L**0.25 * S_f**0.302 * (1+K_ws)**0.04 * (L/D)**0.10

    W_mainLandingGear = 0.0106* K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**(-0.5) * V_stall**0.1

    W_noseLandingGear = 0.032 * K_np * W_l**0.646 * N_l**0.2 * L_n**0.5 * N_nw **0.45

    W_nacelleGroup = 0.6724 * K_ng * N_Lt**0.10 * N_w**0.294 * N_z**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224

    W_engineControls = 5.0*N_en + 0.80*L_ec

    W_starterPneumatic = 49.19*(N_en*W_en/1000)**0.541

    W_fuelSystem = 2.405 * V_t**0.606 * (1 + V_i/V_t)**(-1) * (1 + V_p/V_t) * N_t**0.5

    W_flightControls = 145.9 * N_f**0.554 * (1 + N_m/N_f)**(-1) * S_cs**0.20 * (I_y * 10**(-6))**0.70

    W_APUInstalled = 2.2 * W_APUUninstalled 

    W_instruments = 4.509 * K_r * K_tp * N_c**0.541 * N_en * (L_f + B_w)**0.5

    W_hydraulics = 0.2673 * N_f * (L_f + B_w)**0.937

    W_electrical = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10

    W_avionics = 1.73 * W_uav**0.983

    W_furnishings = 0.0577*N_c**0.1 * W_c**0.393 * S_f**0.75

    W_airConditioning = 62.36*N_p**0.25*(V_pr/1000)**0.604*W_uav**0.10

    W_antiIce = 0.002*W_dg

    W_handling_gear = 0.0003 * W_dg 

    # Miscallaneous weights not accounted by Raymer
    W_cabinFurnishings = 28*(N_p+N_c-3) + 5000 # simplified estimation of cabin furnishings weight
    W_lavatoriesGalleys = 2*300 + 4* 600 # 2 lavatories at 300 lb each, 4 galleys at 600 lb each
    W_minimumFluids = 300 + 600 + 500 # hydraiulic fluid, engine oil, unusable fuel

    W_engines = N_en * W_en
    
    W_misc = 0.03 * W_dg  # paint, insulation, etc.

    
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
        W_engines,
        W_cabinFurnishings,
        W_lavatoriesGalleys,
        W_minimumFluids,
        W_misc
    ]
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
                "Engines",
                "Cabin Furnishings",
                "Lavatories and Galleys",
                "Minimum Fluids",
                "Misc. (paint, insulation, etc.)"
    ]
    for i in range(len(subsystem_names)):
        print(subsystem_names[i], "Weight:", subsystem_values[i], " lb")

    OEM = sum(subsystem_values)
    OEM_SI = OEM * 0.45359237  # Convert to kg
    MTOW_SI = (OEM_SI+18960)*1.459328551 #Convert to MTOW using the fuel mass fraction obtained from the Breguet range equation for jet aircraft
    fuel_mass_SI = MTOW_SI - OEM_SI - 18960 

    #Create a bar chart to visualize the weight breakdown
    #plt.figure(figsize=(12, 6))
    #plt.bar(subsystem_names, subsystem_values)
    #plt.xlabel('Subsystems')
    #plt.ylabel('Weight (lb)')
    #plt.title('Aircraft Weight Breakdown by Subsystem')
    #plt.xticks(rotation=45, ha='right')
    #plt.tight_layout()
    #plt.show()
    
    return subsystem_values


def plotWeightBreakdown(subsystem_values_lb):
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
                "Engines",
                "Cabin Furnishings",
                "Lavatories and Galleys",
                "Minimum Fluids",
                "Misc. (paint, insulation, etc.)"
    ]

    subsystem_values_firstIteration = subsystem_values_lb[0]
    subsystem_values_lastIteration = subsystem_values_lb[-1]  
    #Create a bar chart to visualize the weight breakdown
    plt.figure(figsize=(12, 6))
    plt.bar(subsystem_names, subsystem_values_firstIteration, label='First Iteration, total weight: {:.2f} lb'.format(sum(subsystem_values_firstIteration)), color='blue')
    plt.bar(subsystem_names, subsystem_values_lastIteration, label='Last Iteration, total weight: {:.2f} lb'.format(sum(subsystem_values_lastIteration)), color='orange')
    plt.xlabel('Subsystems')
    plt.ylabel('Weight (lb)')
    plt.title('Aircraft Weight Breakdown by Subsystem')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.legend()
    plt.show()
    
def getClassIMTOW(LiftDragRatio, OEM_kg):
    bypassRatio = 6.0 # PW2040D value
    TSFC = 22*bypassRatio**(-0.19)
    e_f = 4.4*10**1 # J/kg
    eta_jet = prm.V_cruise/(TSFC*e_f)
    print(eta_jet)

    equivalentRange_lst = [11867, 14147, 15177] #[km]
    payload_lst = [18960, 8531, 0]
    MTOW_lst = [] 
    for i in range(3):
        flightMassFraction = np.exp(equivalentRange_lst[i]*1000/(eta_jet*LiftDragRatio*(4.4*10**7/9.81)))
        print("For a range of", equivalentRange_lst[i], "km, the flight mass fraction is:", flightMassFraction)
        MTOW = (OEM_kg+payload_lst[i])*flightMassFraction

        MTOW_lst.append(MTOW)
    MTOW = max(MTOW_lst)
    return MTOW

getClassIMTOW(LiftDragRatio=1, OEM_kg=1)