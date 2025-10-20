import math 
import numpy as np
import params as prm



F_w =       # Fuselage width at horizontal tail intersection, [ft]
Lambda_vt = # horizontal tail sweep at 1/4 (25%) MAC [deg]
t_c_root =      #
lambda_ =       # ???? Taper ratio [-]
Lambda_ht =         # horizontal tail sweep at 1/4 (25%) MAC [deg]
N_l =       # ultimate landing load factor; =N_gear*1.5 [-]


### Constant parameters ###
A_h = prm.empg_AR_h       # ???? Horizontal tail Aspect ratio, [-]
A_v = prm.empg_AR_v     # Vertical Tail Aspect ratio, [-]
K_door = 1.06       # For one side cargo door
K_Lg = 1.0      #1.12 if fuselage-mounted main landing gear; =1.0 otherwise
K_ng = 1.017        # for pylon-mounted nacelle
K_np, K_mp = 1.0, 1.0      #1.126 for kneeling gear; =1.0 otherwise
K_r = 1.0       # 1.133 for reciprocating engine 
K_tp = 1.0      # 0.793 if turboprop
K_uht = 1.0 # 1.143 for unit (all-moving) horizontal tail; 1.0 otherwise
H_t_over_H_v = 0        #0.0 for conventional tail, [-]
N_c =       # number of crew, [-]
N_en =      # number of engines, [-]
N_f = 7     # number of functions performed by controls (typically 4-7)
N_Lt =      # nacelle length, [ft]
N_mss =     # number of main gear shock struts, [-]
N_m = 0     # number of mechanical functions (typically 0-2)
N_mw =      # number of main wheels, [-]
N_nw =      # number of nose wheels, [-]
N_t =       # number of fuel tanks, [-]
N_w =       # nacelle width, [-]
S_n =       # nacelle wetted area, [ft^2]
R_kva = 55  # system electrocal rating, typical values for cargo aircrafts
V_i =       # integral tanks volume, [gal]
V_t =       # total fuel volume, [gal]
V_p = 0     # self-sealing "protected" tanks volume, [gal], apparently only military aircraft
W_en =      # engine weifht, each, [lb]
W_fw =      # weight of fuel in wing, [lb]
W_uav = 1200    #uninstalled avionics weight, [lb]




### Input parameters ###
A=    # aspect ratio [-]
B_w =       # wing span, [ft]
B_h =       # horizontal tail span, [ft]
D =         # fuselage structural depth, [ft]
L =         # fuselage structural length (excludes randome, tail cap), [ft]§
L_f =       # total fuselage length, [ft] ????
L_t =       # tail length, wing quater-MAC to tail-quater-MAC, [ft]
L_m =       # length of main landing gear, [in]
L_n =       # length of nose lanfing gear, [in]
N_z=        # ultimate load factor; 1.5* limit load factor, [-]
S_csw =     # control surface area (wing-mounted), [ft^2]

S_f =       # Fuselage wetted area, [ft^2]
S_ht =      # Horizontal tail area, [ft^2]
S_vt =      # Vertical tail area, [ft^2]
S_w=        #trapezoidal wing area, [ft^2]

V_stall =       #???? Stall speed, [ft/s]?????

W_dg=       # Design gross weight, [lb]
W_l =       # Landing design gross weight, [ln]

Lambda =        # wing sweep at 1/4 (25%) MAC [deg]


### Parameters based on input parameters ###
K_y = 0.3*L_t      #aircraft pitching radius of gyration, [ft] (approx 0.3 L_t)
K_z = 1*L_t        # Aircraft Yawing radius of gyration, [ft] (approx L_t)
S_e =       #elevator area, [ft]
W_ec = 2.331*W_engine**0.901*K_p*K_tr     # weight of engine and contents, (per nacelle), [lb]
K_ws = 0.75*((1+2*lambda_)/(1+lambda_))*(B_w * np.tan(np.radians(Lambda))/L)
N_gen = N_en    # number of generators (typically =N_en)
L_a =       # electrical routing distance, generators to avionics to cockpit, [ft]



### Functions ###
W_wing = 0.0051*(W_dg*N_z)**0.557 * S_w**0.649 * A**0.5 * t_c_root**(-0.4) * (1+lambda_)**0.1 * np.cos(np.radians(Lambda))**(-1)*S_csw**0.1

W_horizontalTail = 0.0379 * K_uht * (1+F_w/B_h)**(-0.25) * W_dg**0.639 * N_z**0.1 * S_ht**0.75 *L_t**(-1) * K_y**0.704 * (np.cos(np.radians(Lambda_ht)))**(-1) * A_h**0.166 * (1+S_e/S_ht)**0.1

W_verticalTail = 0.0026 * (1+H_t_over_H_v)**0.225 * W_dg**0.556 * N_z**536 * L_t**(-0.5) * S_vt**0.5 * K_z**0.875 * np.cos(np.radians(Lambda_vt))**(-1) * A_v**0.35 * (t_c_root)**(-0.5)

W_fuselage = 0.3280 * K_door * K_Lg *(W_dg*N_z)**0.5 * L**0.25 * S_f**0.302 * (1+K_ws)**0.04 * (L/D)**0.10

W_mainLandingGear = 0.0106* K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**(-0.5) * V_stall**0.1

W_noseLandingGear = 0.032 * K_np * W_l**0.646 * N_l**0.2 * L_n**0.5 * N_nw **0.45

#W_installedEngineTotal = 2.575 * W_en**0.922 * N_en

W_nacelleGroup = 0.6724 * K_ng * N_Lt**0.10 * N_w**0.294 * N_z**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224

W_engineControls = 5.0*N_en + 0.80*L_ec

W_starterPneumatic = 49.19*(N_en*W_en/1000)**0.541

W_fuelSystem = 2.405 * V_t**0.606 * (1+V_i/V_t)*(-1) * (1 + V_p/V_t)*N_t**0.5

W_flightControls = 145.9 * N_f**0.554 * (1 + N_m/N_f)**(-1) * S_cs**0.20 * (I_y * 10**(-6))*0.07¨

W_APUInstalled = 2.2 * W_APUUninstalled 

W_instruments = 4.509 * K_r * K_tp * N_c**0.541 * N_en * (L_f + B_w)**0.5

W_hydraulics = 0.2673 * 0.2673 * N_f * (L_f + B_w)**0.937

W_electrical = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10

W_avionics = 1.73 * W_uav**0.983

