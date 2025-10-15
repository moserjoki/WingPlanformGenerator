import math
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# ---------------------------
# Airfoil Class
# ---------------------------
class Airfoil:
    def __init__(self, filepath, c_l_alpha, c_d0, tau, c_l_max):
        self.filepath = filepath
        self.x, self.y, self.x_upper, self.y_upper, self.x_lower, self.y_lower = self._import_airfoil(filepath)
        self.c_l_alpha = c_l_alpha
        self.c_d0 = c_d0
        self.tau = tau
        self.c_l_max = c_l_max

    def _import_airfoil(self, filepath):
        x, y = [], []
        with open(filepath, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) == 2:
                    try:
                        x_val, y_val = map(float, parts)
                        x.append(x_val)
                        y.append(y_val)
                    except ValueError:
                        continue

        x = np.array(x)
        y = np.array(y)

        # --- Find leading edge (minimum x) ---
        le_index = np.argmin(x)

        # Split into upper and lower surfaces
        x_upper = x[:le_index + 1]
        y_upper = y[:le_index + 1]
        x_lower = x[le_index:]
        y_lower = y[le_index:]

        # Ensure both are ordered left-to-right (increasing x)
        if x_upper[0] > x_upper[-1]:
            x_upper = np.flip(x_upper)
            y_upper = np.flip(y_upper)
        if x_lower[0] > x_lower[-1]:
            x_lower = np.flip(x_lower)
            y_lower = np.flip(y_lower)

        return x, y, x_upper, y_upper, x_lower, y_lower
    
    def get_thickness_at(self, x_c):
        y_u = np.interp(x_c, self.x_upper, self.y_upper)
        y_l = np.interp(x_c, self.x_lower, self.y_lower)
        return y_u - y_l
    
    def cross_section(self, xc_start, xc_end, dxc, plot):
        """
        Numerically integrates the airfoil thickness distribution between
        xc_start and xc_end (chordwise positions) to calculate spar positions.
        Also stores thickness values for visualization.
        """
        xc = xc_start
        sum_val = 0.0

        while xc < xc_end:
            t = self.get_thickness_at(xc)
            sum_val += t*dxc

            if plot:
                # Interpolate along the stored lower surface
                y_lower = np.interp(xc, self.x_lower, self.y_lower)
                y_upper = np.interp(xc, self.x_upper, self.y_upper)
                plt.fill_between([xc, xc + dxc], y_lower, y_upper, color='red', alpha=0.3)

            xc += dxc

        if plot:
            plt.xlabel("x/c []")
            plt.ylabel("y/c []")
            plt.title("Airfoil cross section []")
            plt.plot(self.x, self.y, 'k')  # plot the airfoil outline
            plt.axis("equal")
            plt.show()

        return sum_val

# ---------------------------
# Matching Diagram Class
# ---------------------------
class MatchingDiagram:
    g = 9.81 # [m/s^2] gravitational acceleration

    WS_array = np.array([]) # [N/m^2]
    TW_cruise_speed_array = np.array([])
    TW_cruise_speed = 0 
    TW_climb_array = np.array([])
    TW_climb = 0 
    TW_cg119_array = np.array([])
    TW_cg119 = 0 
    TW_cg121A_array = np.array([])
    TW_cg121A = 0 
    TW_cg121B_array = np.array([])
    TW_cg121B = 0 
    TW_cg121C_array = np.array([])
    TW_cg121C = 0 
    TW_cg121D_array = np.array([])
    TW_cg121D = 0 
    TW_TOFL_array = np.array([])

    WS_min_speed = 0
    WS_landing_field = 0

    def __init__(self, m_MTOW, AR, V_cruise, M_cruise, ρ_cruise):
        self.m_MTOW = m_MTOW
        self.AR = AR
        self.V_cruise = V_cruise
        self.M_cruise = M_cruise
        self.ρ_cruise = ρ_cruise

        self.TW_min = None
        self.WS_min = None

    def compute(self, C_L_max_takeoff, C_L_max_landing):

        self.WS_array = np.array([]) # [N/m^2]
        self.TW_cruise_speed_array = np.array([])
        self.TW_cruise_speed = 0 
        self.TW_climb_array = np.array([])
        self.TW_climb = 0 
        self.TW_cg119_array = np.array([])
        self.TW_cg119 = 0 
        self.TW_cg121A_array = np.array([])
        self.TW_cg121A = 0 
        self.TW_cg121B_array = np.array([])
        self.TW_cg121B = 0 
        self.TW_cg121C_array = np.array([])
        self.TW_cg121C = 0 
        self.TW_cg121D_array = np.array([])
        self.TW_cg121D = 0 
        self.TW_TOFL_array = np.array([])

        self.WS_min_speed = 0
        self.WS_landing_field = 0

        V_app = 73 # [m/s] Approach speed (chosen)
        beta = 0.82 # [] Mass fraction 

        l_landing_field = 1981.2 # [m] length of landing field for landing req
        h_landing = 0 # [m] altitude for landing req
        ΔT_landing = 0 # [K] temperature difference for landing req
        ρ_landing = 1.225 # [kg/m^3] Density 
        C_lfl = 0.45 # [] Landing field length coefficient

        self.WS_min_speed = (1/beta)*(ρ_landing/2)*(V_app/1.23)**2*C_L_max_landing # [N/m^2]
        self.WS_landing_field = (1/beta)*ρ_landing*(l_landing_field/C_lfl)*(C_L_max_landing/2) # [N/m^2]

        # CRUISE SPEED REQUIREMENT 
        α_P = 0.228 # [] Power/thrust lapse 
        beta_cr = 0.95 # [] mass fraction for the cruise speed requirement 
        c_D0 = 0.017 # [] zero lift drag coefficient 
        e = 0.74 # [] oswald efficiency factor

        # CLIMB RATE REQUIREMENT 
        ρ_climb = 0.56 # [kg/m^3]
        T_climb = 240 # [K]
        p_climb = 38800 # [Pa] 
        p_SL_ISA = 101325 # [Pa] sea level pressure
        T_SL_ISA = 288.15 # [K]
        R = 287 # [J/(kg*K)] ideal gas constant air 
        γ = 1.4 # [] ratio of specific heats 
        B = 8 # [] bypass ratio
        c = 6.2 # [m/s] climb rate 
        beta_climb = 0.95 # [] mass fraction 

        # CLIMB GRADIENT 119
        beta_cg119 = 1 # [] 
        ΔT_cg119 = 0 # [K]
        coverV_cg119 = 3.2 # [%] climb gradient c/V 
        c_D0_cg119 = 0.084
        e_cg119 = 0.86
        T_cg119= 288.15  # [K]
        ρ_cg119 = 1.225  # [kg/m^3]

        # CLIMB GRADIENT 121A
        beta_cg121A = 1 # [] 
        ΔT_cg121A = 0 # [K]
        coverV_cg121A = 0 # [%] climb gradient c/V 
        c_D0_cg121A = 0.058
        e_cg121A = 0.81
        T_cg121A = 288.15  # [K]
        ρ_cg121A = 1.225 # [kg/m^3]

        # CLIMB GRADIENT 121B
        beta_cg121B = 1 # [] 
        ΔT_cg121B = 0 # [K]
        coverV_cg121B = 2.4 # [%] climb gradient c/V 
        c_D0_cg121B = 0.039
        e_cg121B = 0.81
        T_cg121B = 288.15  # [K]
        ρ_cg121B  = 1.225 # [kg/m^3]

        # CLIMB GRADIENT 121C
        beta_cg121C = 1 # []
        ΔT_cg121C = 0 # [K]
        coverV_cg121C = 1.2 # [%] climb gradient c/V 
        c_D0_cg121C = 0.019
        e_cg121C = 0.77
        T_cg121C = 288.15  # [K]
        ρ_cg121C = 1.225 # [kg/m^3]

        # CLIMB GRADIENT 121D
        beta_cg121D = 0.84 # []
        ΔT_cg121D = 0 # [K]
        coverV_cg121D = 2.1 # [%] climb gradient c/V 
        c_D0_cg121D = 0.065
        e_cg121D = 0.86
        T_cg121D = 288.15  # [K]
        ρ_cg121D = 1.225 # [kg/m^3]

        # TAKE-OFF FIELD LENGTH - ALL ENGINES OPERATIVE
        L_TOFL = 3048 # [m]
        h_TOFL = 500 # [m]
        ΔT_TOFL = 15 # [K]
        c_D0_TOFL = 0.017 # [K]
        e_TOFL = 0.74 # []
        T_TOFL = 298 # [K]
        ρ_TOFL = 1.11 # !A [kg/m^3]
        p_TOFL = 95500	# !A [Pa]
        k_T_TOFL = 0.85 # []
        h_2_TOFL = 11 # [m]
        #print(T_TOFL)

        WS = 100
        WS_step = 100 
        WS_max = max(self.WS_min_speed, self.WS_landing_field) + 600 # Calculates the diagram limits by taking the furthest right wing loading

        TW_TOFL = 0 

        while WS < WS_max:
            TW_cruise_speed = ((beta_cr/α_P)*((c_D0*(1/2)*self.ρ_cruise*self.V_cruise**2)/(beta_cr*WS) + (beta_cr*WS)/(math.pi*self.AR*e*(1/2)*self.ρ_cruise*self.V_cruise**2)))
            self.TW_cruise_speed_array = np.append(self.TW_cruise_speed_array, TW_cruise_speed)

            # CRUISE SPEED REQUIREMENT 
            M_climb = math.sqrt(((WS*2)/(ρ_climb*C_L_max_takeoff))/(R*T_climb*γ))
            p_t_climb = p_climb*(1+0.2*M_climb**2)**(1.4/0.4)
            δ_t_climb = p_t_climb/p_SL_ISA
            α_T_climb = δ_t_climb*(1-(0.43+0.014*B)*math.sqrt(M_climb))
            TW_climb = (beta_cr/α_T_climb)*(math.sqrt(((c**2)/(beta_climb*WS))*(ρ_climb/2)*math.sqrt(c_D0*math.pi*self.AR*e)) + 2*math.sqrt(c_D0/(math.pi*self.AR*e)))
            self.TW_climb_array = np.append(self.TW_climb_array, TW_climb)

            # CLIMB GRADIENT 119
            M_cg119 = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_cg119*R*T_cg119*γ)))
            p_t_cg119 = p_SL_ISA*(1+0.2*M_cg119**2)**(1.4/0.4)
            δ_cg119 = p_t_cg119/p_SL_ISA
            α_T_cg119 = δ_cg119*(1-(0.43+0.014*B)*math.sqrt(M_cg119))
            TW_cg119 = (beta_cg119/α_T_cg119)*(coverV_cg119*0.01 + 2*math.sqrt(c_D0_cg119/(math.pi*self.AR*e_cg119)))
            self.TW_cg119_array = np.append(self.TW_cg119_array, TW_cg119)

            # CLIMB GRADIENT 121A
            M_cg121A = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_cg121A*R*T_cg121A*γ)))
            p_t_cg121A = p_SL_ISA*(1+0.2*M_cg121A**2)**(1.4/0.4)
            δ_cg121A = p_t_cg121A/p_SL_ISA
            α_T_cg121A = δ_cg121A*(1-(0.43+0.014*B)*math.sqrt(M_cg121A))
            TW_cg121A = 2*(beta_cg121A/α_T_cg121A)*(coverV_cg121A*0.01 + 2*math.sqrt(c_D0_cg121A/(math.pi*self.AR*e_cg121A)))
            self.TW_cg121A_array = np.append(self.TW_cg121A_array, TW_cg121A)

            # CLIMB GRADIENT 121B
            M_cg121B = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_cg121B*R*T_cg121B*γ)))
            p_t_cg121B = p_SL_ISA*(1+0.2*M_cg121B**2)**(1.4/0.4)
            δ_cg121B = p_t_cg121B/p_SL_ISA
            α_T_cg121B = δ_cg121B*(1-(0.43+0.014*B)*math.sqrt(M_cg121B))
            TW_cg121B = 2*(beta_cg121B/α_T_cg121B)*(coverV_cg121B*0.01 + 2*math.sqrt(c_D0_cg121B/(math.pi*self.AR*e_cg121B)))
            self.TW_cg121B_array = np.append(self.TW_cg121B_array, TW_cg121B)
            
            # CLIMB GRADIENT 121C
            M_cg121C = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_cg121C*R*T_cg121C*γ)))
            p_t_cg121C = p_SL_ISA*(1+0.2*M_cg121C**2)**(1.4/0.4)
            δ_cg121C = p_t_cg121C/p_SL_ISA
            α_T_cg121C = δ_cg121C*(1-(0.43+0.014*B)*math.sqrt(M_cg121C))
            TW_cg121C = 2*(beta_cg121C/α_T_cg121C)*(coverV_cg121C*0.01 + 2*math.sqrt(c_D0_cg121C/(math.pi*self.AR*e_cg121C)))
            self.TW_cg121C_array = np.append(self.TW_cg121C_array, TW_cg121C)

            # CLIMB GRADIENT 121D
            M_cg121D = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_cg121D*R*T_cg121D*γ)))
            p_t_cg121D = p_SL_ISA*(1+0.2*M_cg121D**2)**(1.4/0.4)
            δ_cg121D = p_t_cg121D/p_SL_ISA
            α_T_cg121D = δ_cg121D*(1-(0.43+0.014*B)*math.sqrt(M_cg121D))
            TW_cg121D = 2*(beta_cg121D/α_T_cg121D)*(coverV_cg121D*0.01 + 2*math.sqrt(c_D0_cg121D/(math.pi*self.AR*e_cg121D)))
            self.TW_cg121D_array = np.append(self.TW_cg121D_array, TW_cg121D)

            # TAKE-OFF FIELD LENGTH - ALL ENGINES OPERATIVE
            M_TOFL = math.sqrt(((WS*2)/(C_L_max_takeoff*ρ_TOFL*R*T_TOFL*γ)))
            p_t_TOFL = p_TOFL*(1+0.2*M_TOFL**2)**(1.4/0.4)
            δ_TOFL = p_t_TOFL/p_SL_ISA
            α_T_TOFL = δ_TOFL*(1-(0.43+0.014*B)*math.sqrt(M_TOFL))
            TW_TOFL = (1/α_T_TOFL)*(1.15*math.sqrt(WS/(L_TOFL*k_T_TOFL*ρ_TOFL*self.g*math.pi*self.AR*e_TOFL))+4*h_2_TOFL/L_TOFL)
            self.TW_TOFL_array = np.append(self.TW_TOFL_array, TW_TOFL)

            self.WS_array  = np.append(self.WS_array, WS) 
            WS = WS + WS_step 

        # FINDING THE DESIGN POINT 
        self.WS_min = min(self.WS_min_speed, self.WS_landing_field)
        WS_search = round(self.WS_min, -2)
        i_search = int(WS_search/100 )

        self.TW_min = max(
            self.TW_climb_array[i_search],
            self.TW_cg119_array[i_search],
            self.TW_cg121A_array[i_search],
            self.TW_cg121B_array[i_search],
            self.TW_cg121C_array[i_search],
            self.TW_cg121D_array[i_search],
            self.TW_TOFL_array[i_search]
        )


        #OUTPUT
        S_w = self.m_MTOW*self.g/self.WS_min

        print("WTO/SW:", round(self.WS_min, 4), "N/m^2")
        print("TWO/WTO:", round(self.TW_min, 4), "N/N")
        print("Wing area:", round(S_w, 4), "m^2")

        return S_w

    def plot(self):
        # optional plotting of matching diagram

        # PLOTTING RESULTS 
        plt.axvline(x=self.WS_min_speed, color="blue", linestyle="-", label="Minimum speed")
        plt.axvline(x=self.WS_landing_field, color="darkorange", linestyle="-", label="Landing field length")
        plt.plot(self.WS_array, self.TW_cruise_speed_array, color="purple", linestyle="-", label="Cruise speed")
        plt.plot(self.WS_array, self.TW_climb_array, color="gold", linestyle="-", label="Climb rate")
        plt.plot(self.WS_array, self.TW_cg119_array, color="crimson", linestyle="-", label="Climb gradient 119")
        plt.plot(self.WS_array, self.TW_cg121A_array, color="red", linestyle="-", label="Climb gradient 121A")
        plt.plot(self.WS_array, self.TW_cg121B_array, color="orange", linestyle="-", label="Climb gradient 121B")
        plt.plot(self.WS_array, self.TW_cg121C_array, color="goldenrod", linestyle="-", label="Climb gradient 121C")
        plt.plot(self.WS_array, self.TW_cg121D_array, color="skyblue", linestyle="-", label="Climb gradient 121D")
        plt.plot(self.WS_array, self.TW_TOFL_array, color="limegreen", linestyle="-", label="Take off field length")

        # FIND DESIGN POINT
        plt.plot(self.WS_min, self.TW_min, marker="x", color="darkgreen", markersize=10, markeredgewidth=2, label="Design point")

        plt.ylim(0, 0.5)  # y-axis from 0 to 30
        plt.xlabel("W/S [N/m^2]")
        plt.ylabel("T/W [N/N]")
        plt.title("Matching Diagram")
        plt.legend()
        plt.grid(True)
        plt.show()
        pass

# ---------------------------
# Wing Sizing Class
# ---------------------------
class WingSizing:
    g = 9.81 # [m/s^2] gravitational acceleration

    def __init__(self, m_MTOW, AR, V_cruise, M_cruise, ρ_cruise, airfoil: Airfoil):
        # Inputs
        self.m_MTOW = m_MTOW
        self.AR = AR
        self.V_cruise = V_cruise
        self.M_cruise = M_cruise
        self.ρ_cruise = ρ_cruise

        # Gather data from airfoil 
        self.c_l_alpha = airfoil.c_l_alpha
        self.c_d0 = airfoil.c_d0
        self.tau = airfoil.tau
        self.c_l_max = airfoil.c_l_max

        # Wing planform parameters
        self.S_w = None
        self.b = None
        self.X_LE = None
        self.X_TE = None
        self.c_root = None
        self.c_tip = None
        self.taper_ratio = None
        self.MAC = None
        self.leading_sweep = None
        self.quarter_sweep = None
        self.trailing_sweep = None
        self.dihedral = None
        self.b1 = None
        self.b2 = None


    # ---------------------------
    # Planform sizing
    # ---------------------------
    def planform_sizing(self, S_w):
        self.S_w = S_w 
        self.b = math.sqrt(self.AR*self.S_w) # [m] wing span

        # Given values
        k_a = 0.935             # [] technology factor for all super critical airfoils
        t_c_streamwise = 0.14   # [] t/c 
        M_DD = self.M_cruise + 0.02 # [] target drag divergence Mach number
        
        C_L_des = 2*1.1*self.m_MTOW*self.g/(self.S_w*self.ρ_cruise*self.V_cruise**2) # []

        quart_sweep_prelim = math.degrees(math.acos(1.16/(self.M_cruise+0.5)))

        M_DD_prelim = k_a/math.cos(math.radians(quart_sweep_prelim))-t_c_streamwise/math.cos(math.radians(quart_sweep_prelim))**2-C_L_des/(10*math.cos(math.radians(quart_sweep_prelim))**3)
       
        # Equation to solve:
        # M_DD = k_a/cos(x) - t_c/cos(x)^2 - C_L/(10*cos(x)^3)
        def equation(x):
            return k_a / np.cos(x) - t_c_streamwise*np.cos(x) / np.cos(x)**2 - C_L_des / (10 * np.cos(x)**3) - M_DD

        # Initial guess (in radians)
        x0 = 0.3   # ~17 degrees

        # Solve
        solution = fsolve(equation, x0)
        x_rad = solution[0]
        x_deg = math.degrees(x_rad)
        self.quart_sweep = x_deg

        self.taper_ratio = 0.2*(2-self.quart_sweep*math.pi/180)
        self.c_root = 2*self.S_w/((self.taper_ratio+1)*self.b)
        self.c_tip = self.c_root*self.taper_ratio
        
        self.leading_sweep = math.degrees(math.atan(math.tan(math.radians(self.quart_sweep))-(self.c_root/(2*self.b))*(self.taper_ratio-1)))
        self.trailing_sweep = math.degrees(math.atan(math.tan(math.radians(self.quart_sweep))+3*(self.c_root/(2*self.b))*(self.taper_ratio-1))) # [deg]
        self.dihedral = 3-0.1*self.quart_sweep+2

        self.MAC = (2/3)*self.c_root*((1+self.taper_ratio+self.taper_ratio**2)/(1+self.taper_ratio))
        self.Y_mac = self.b/6*(1+2*self.taper_ratio)/(1+self.taper_ratio)
        self.X_mac = self.Y_mac*math.tan(math.radians(self.leading_sweep))

        AR_bound = 17.7*(2 - self.taper_ratio)*math.e**(-0.043*self.quart_sweep)

        print(f"AR {self.AR:0.3f} <= AR_bound: {AR_bound:0.3f}")
        print(f"b: {round(self.b, 4)} m | Λ_c/4: {round(self.quart_sweep,4)} ° | Λ_LE: {round(self.leading_sweep,4)} ° | Λ: {round(self.taper_ratio,4)} |  c_r: {round(self.c_root,4)} m | c_t: {round(self.c_tip,4)} m | dihedral: {self.dihedral:0.2f} °")
        print(f"MAC: {round(self.MAC,4)} m | X_mac: {self.X_mac:0.3f} m | Y_mac: {self.Y_mac:0.3f} m")
        return

    # ---------------------------
    # Helper function
    # ---------------------------
    def leading_edge(self, y):
        return math.tan(math.radians(self.leading_sweep)) * y

    def trailing_edge(self, y):
        return math.tan(math.radians(self.trailing_sweep)) * y + self.c_root

    def chord(self, y):
        return self.trailing_edge(y) - self.leading_edge(y)

    # ---------------------------
    # Plotting method
    # ---------------------------
    def plot(self):
        y = np.linspace(0, self.b/2, 200)
        LE = np.array([self.leading_edge(yi) for yi in y])
        TE = np.array([self.trailing_edge(yi) for yi in y])

        # Fill area first
        plt.fill_between(y, -LE, -TE, color='gray', alpha=0.3)

        # Plot lines on top
        plt.plot(y, -LE, label='Leading Edge', color='gray', linewidth=2)
        plt.plot(y, -TE, label='Trailing Edge', color='gray', linewidth=2)

        plt.xlabel('y [m]')
        plt.ylabel('x [m]')
        plt.title('Wing Planform')
        plt.legend()

        # Fixed axes
        plt.axis('equal')
        plt.xlim(0, max(y)*1.1)
        plt.ylim(-max(TE)*1.1,0)

        plt.show()

    # ---------------------------
    # Aileron sizing
    # ---------------------------
    def aileron_sizing(self, C_L_max_landing):
        delta_a_up = 17
        f_differential_aileron = 0.75
        delta_a_down = delta_a_up*f_differential_aileron
        delta_a = (1/2)*(delta_a_up+delta_a_down)

        ρ_SL = 1.22522568 # [kg/m^3]
        
        V_stall = math.sqrt(2*self.m_MTOW*self.g/(ρ_SL*C_L_max_landing*self.S_w)) # [m/s]
        print(f"V_stall: {V_stall:0.3f}")
        delta_bank_req = 60  # deg
        delta_t_req = 7 # s

        b2_percent = 0.75            # Most outward position of ailerons (not tip to minimize alerion reversal)
        self.b2 = b2_percent * self.b / 2

        # ---------------------------
        # FUNCTIONS
        # ---------------------------
        def integral_aileron_control(b1, b2):
            sum_val = 0
            dy = 0.001
            y = b1
            while y < b2:
                sum_val += y*self.chord(y)*dy
                y += dy
            return sum_val

        def integral_roll_damping(b1, b2):
            sum_val = 0
            dy = 0.001
            y = b1
            while y < b2:
                sum_val += y**2*self.chord(y)*dy
                y += dy
            return sum_val

        def delta_t_from_b1(b1, printing):
            C_l_delta_alpha = ((2*self.c_l_alpha*self.tau)/(self.S_w*self.b)) * integral_aileron_control(b1, self.b2)
            C_l_p = (-(4*(self.c_l_alpha+self.c_d0))/(self.S_w*self.b**2)) * integral_roll_damping(0, self.b/2)
            P = -(C_l_delta_alpha / C_l_p) * math.radians(delta_a) * (2*V_stall / self.b)
            delta_t = math.radians(delta_bank_req) / P

            if printing:
                print(f"C_l_delta_alpha: {C_l_delta_alpha :0.3f} | C_l_p: {C_l_p:0.3f} ")
                print(f"P: {P:0.3f} rad/s | {math.degrees(P):0.3f} °/s")
                print(f"delta_t: {delta_t:0.3f} s (<= {delta_t_req} s) to turn delta_bank_req: {delta_bank_req} °")
            return delta_t - delta_t_req  # we want this = 0

        # ---------------------------
        # SOLVE FOR b1
        # ---------------------------
        sol = root_scalar(lambda b1: delta_t_from_b1(b1, printing=False), bracket=[0.01, self.b2-0.01], method='bisect')

        if sol.converged:
            self.b1 = sol.root
        
        delta_t_from_b1(self.b1, True)
        print(f"b2: {self.b2:0.4f} m | {b2_percent:0.4f}% of b/2 | b1: {self.b1:0.4f} m | {self.b1/(self.b/2):0.4f}% of b/2 | b2-b1: {(self.b2-self.b1):0.2f} m")
        return
    
    # ---------------------------
    # Transforms airfoil c_l_max to wing C_L_max for clean configuration
    # ---------------------------
    def DATCOM_C_L_max_clean(airfoil: Airfoil, C_L_C_l_max_ratio, delta_C_L_max):
        C_L_max_airfoil = C_L_C_l_max_ratio*airfoil.c_l_max
        C_L_max_clean =  C_L_max_airfoil + delta_C_L_max
        return C_L_max_clean

    # ---------------------------
    # High-Lift Devices sizing
    # ---------------------------
    def HLD_sizing(self, C_L_max_clean):    
        # Wing planform details: 
        d_fuselage = 5.2

        # Leading Edge Fraction Flapped Area:
        S_wfLE_to_S = 0.8

        # Find the reference flapped area for TE devices
        S_wfTE_half = (self.c_root*self.b1-(self.c_root-self.c_tip)/(self.b)*(self.b1)**2)-(self.c_root*(d_fuselage/2)-(self.c_root-self.c_tip)/(self.b)*((d_fuselage/2))**2)
        S_wfTE_to_S = (2/(self.S_w))*S_wfTE_half

        print(f"S_wfLE_to_S: {S_wfLE_to_S:0.2f} | S_wfTE_to_S: {S_wfTE_to_S:0.2f}")

        # Lift coefficient of reference clean wing:
        # Lift contribution of leading edge devices to delta_C_l:
        delta_C_l_LE = 0.4
        # Lift contribution of trailing edge devices to delta_C_l, as tuple (Take-off, Landing):
        delta_C_l_TE = (1.88, 1.94)

        ### Find contribution of LE devices to delta_C_L
        delta_C_L_LE = 0.9 * delta_C_l_LE * S_wfLE_to_S * math.cos(math.radians(self.leading_sweep))

        ### Find contributuion of TE devices to delta_C_L for (take-off and landing)
        delta_C_L_TE = (0.9 * delta_C_l_TE[0] * S_wfTE_to_S * math.cos(math.radians(self.trailing_sweep)), 0.9 * delta_C_l_TE[1] * S_wfTE_to_S * math.cos(math.radians(self.trailing_sweep)))

        delta_C_L = (C_L_max_clean  + delta_C_L_LE + delta_C_L_TE[0], C_L_max_clean + delta_C_L_LE + delta_C_L_TE[1])
        print(f"C_L_max_clean: {C_L_max_clean:0.2f} | C_L_max_takeoff: {delta_C_L[0]:0.4f} | C_L_max_landing: {delta_C_L[1]:0.4f}")

        return delta_C_L[0], delta_C_L[1]
    

    # ---------------------------
    # Calculate Wing Volume
    # ---------------------------
    def volume_section(self, b_start, b_end, S_ref_chord, db):
        bcur = b_start
        volume = 0.0
        while bcur < b_end:
            volume += S_ref_chord*self.chord(bcur)**2*db
            bcur += db
        return volume
    
    def empenage_sizing(self, l_fus, X_cg_aft, printing):
        #SELECTED CONSTANTS (all values were selected to be in the middle of
        # acceptable range:
        Vh = 1.01
        Vv = 0.079
        AR_v = 1.5
        AR_h = 4
        taper_v = 0.5
        taper_h = 0.65

        #VERTICAL TAIL:
        #Leading edge sweep is equal to the wing leading edge sweep
        Leading_Edge_Sweep_V=37.3 #COULD CHANGE
        lv = 0.9*l_fus-X_cg_aft
        S_v = (Vv*self.S_w*self.b)/lv
        b_v = math.sqrt(AR_v*S_v)
        c_r_v = 2*S_v/((1+taper_v)*b_v)
        c_t_v = c_r_v*taper_v
        MAC_v = 2/3*c_r_v*((1+taper_v+taper_v**2)/(1+taper_v))

        #HORIZONTAL TAIL:
        #equal to the wing quarter chord or higher (till 40 deg)
        Quarter_Chord_Sweep_H = 34.21 #COULD CHANGE
        lh = 0.9*l_fus-X_cg_aft
        S_h = Vh*self.S_w*self.MAC/lh
        b_h = math.sqrt(AR_h*S_h)
        c_r_h = 2*S_h/((1+taper_h)*b_h)
        c_t_h = c_r_h*taper_h
        MAC_h = 2/3*c_r_h*((1+taper_h+taper_h**2)/(1+taper_h))

        if printing:
            print("VERTICAL TAIL:")
            print("Moment arm:", round(lv,2))
            print("Vertical Tail Area:", round(S_v,2))
            print("Vertical Tail Span:", round(b_v,2))
            print("Root Chord:", round(c_r_v,2))
            print("Tip Chord:", round(c_t_v, 2))
            print("MAC Vertical Tail:", round(MAC_v,2))
            print("")
            print("HORIZONTAL TAIL:")
            print("Moment arm:", round(lh,2))
            print("Horizontal Tail Area:",round(S_h,2))
            print("Horizontal Tail Span:", round(b_h,2))
            print("Root Chord:", round(c_r_h,2))
            print("Tip Chord:", round(c_t_h,2))
            print("MAC Horizontal Tail:", round(MAC_h,2))
        return