import math
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from sizing import *


m_MTOW = 114960 # [kg] Maximum Take Off Weight from Class I Weight Estimation . 
AR = 6.9 # [] aspect ratio

# Create objects used for subsequent calculations
airfoil = Airfoil("airfoils/NASA SC(2)-0414.dat") 
wing = WingSizing(m_MTOW, AR, M_cruise)
cruise_matching_diagram = MatchingDiagram(m_MTOW, AR, M_cruise)


# Dummy values for first matching diagram calculation
C_L_max_take_off_cur = 2.1 # [] Cl during take-off
C_L_max_landing_cur = 2.4 # [] Cl during landing

# Iterating on the design based on the initial given configuration until it converges
for i in range(5):
    print(f"Iteration {i+1}")
    S_w_cur  = cruise_matching_diagram.compute(C_L_max_take_off_cur, C_L_max_landing_cur)
    wing.planform_sizing(S_w_cur)
    wing.aileron_sizing(C_L_max_landing_cur)

    C_L_max_clean = wing.DATCOM_C_L_max_clean()
    C_L_max_take_off_cur, C_L_max_landing_cur = wing.HLD_sizing(C_L_max_clean)
 
cruise_matching_diagram.plot()

X_cg_aft = 21.98 #RANDOM INITIAL VALUE
wing.empenage_sizing(X_cg_aft, True)

wing.fuel_volume(airfoil)

wing.plot()