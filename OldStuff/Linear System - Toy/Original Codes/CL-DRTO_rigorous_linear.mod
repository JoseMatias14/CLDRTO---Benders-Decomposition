
param n_DRTO integer >0;		   	# prediction horizon of DRTO
set ns ordered;                     # Number of scenarios
param prob_scenario{ns} >= 0; 
check 0.99999 <= sum{i in ns} prob_scenario[i] <= 1.00001 ;            # Probability of a given scenario
#--------------------------------------------------------------------------------
# Model parameters
#-------------------------------------------------------------------------------
param nx_;                 	        # No. of differential states
param ny_;                 	        # No. of outputs
param nu_;                 	        # No. of inputs
param A_{1..nx_, 1..nx_};  	        # A matrix for MPC
param B_{1..nx_, 1..nu_};  	        # B matrix for MPC
param C_{1..ny_, 1..nx_};  	        # C matrix for MPC
param A_DRTO{ns, 1..nx_, 1..nx_};  	# A matrix for DRTO
param B_DRTO{ns, 1..nx_, 1..nu_};  	# B matrix for DRTO
param C_DRTO{ns, 1..ny_, 1..nx_};  	# C matrix for DRTO
set grades ordered;                 # Number of product grades
param uss{grades};                  # Control value at steady state operation for each grade
param flow_rate;

param dt_MPC;                                 # Sampling time MPC
param dt_DRTO;                                # Sampling time 
param n_DRTOn = n_DRTO;           # n in the rigorous approach
param sph >0 integer;                         # Control/output values are kept const per N control moves

set N_DRTO  := 0..n_DRTO-1;        # RTO prediction
set N_DRTOn := 0..n_DRTOn-1;       # Prediction horizon for first level of MPC subproblems 
set s_nx    := 1..nx_ ;            # state variables
set s_ny    := 1..ny_;             # output variables
set s_nu    := 1..nu_;             # control variables

param lin_point_yrto{s_ny};        # Linearization point for yrto
param lin_point_urto{s_nu};
param lin_point_xrto{s_nx};
param theta;                                           # Solver feasibility tolerance

param ymeas{s_ny};                    # Measurement for DRTO
param xrto_0{ns,s_nx};                   # Initial value of the state x
param urto_minus_one{s_nu};           # Value of control at time minus one
#----------------------------------------------------------------------------------
# Bounds
#----------------------------------------------------------------------------------
param umin{s_nu};        # Lower bound on control values
param umax{s_nu};        # Upper bound on control values
param ymin{s_ny};        # Lower bound on output values
param ymax{s_ny};        # Upper bound on output values
param ymin_ref{s_ny};    # Lower bound on reference trajectory
param ymax_ref{s_ny};    # Upper bound on reference trajectory
param umin_ref{s_nu};    # Lower bound on reference trajectory
param umax_ref{s_nu};    # Upper bound on reference trajectory

#----------------------------------------------------------------------------------
# DRTO variables
#---------------------------------------------------------------------------------
var urto{ns,s_nu,0..n_DRTO-1};                                                      # Control variables for DRTO
var xrto{ns,s_nx,0..n_DRTO};                                                        # State variables for DRTO
var yrto{ns,s_ny,1..n_DRTO};                                                        # Output variables for DRTO
param d_drto{s in ns, i in s_ny} = ymeas[i] - sum{n in s_nx} C_[i,n]*xrto_0[s, n];  # Constant disturbance               
                            
#******************************************************************************************
# MPC subproblems
#******************************************************************************************
param p_mpc;                              # Prediction horizon in MPC
param m_mpc;                              # Control horizon in MPC
param Q_{1..ny_};                         # Weight matrix output
param R_{1..nu_};                         # Weight matrix move suppresion
param S_{1..nu_};                         # Weight matrix controls
param xmpc_0{i in 1..nx_};                # x0 for MPC

set P_mpc := 0..p_mpc-1;
set M_mpc := 0..m_mpc-1;

# Reference Trajectories
var yref{ns,s_ny, j in 1..n_DRTO + p_mpc - 1};      # Reference trajectory
var uref{ns,s_nu, j in 0..n_DRTO -1 + m_mpc -1};    # Reference trajectory
var bref{ns,s_nu, j in 1..n_DRTO} binary;
param H_ref integer > 0;                                # Number of time-steps the reference trajectory should staty constant
param N_ref integer > 0;                            # Allowed number of changes in ref traj

param M_comp{i in s_nu} default 20*(umax[i]-umin[i]) >=0;                             # For Big M formulation for complementarity constraints

# Variables MPC subproblems
var xmpc{ns, N_DRTOn,s_nx,0..p_mpc};                 # State variables mpc level            
var umpc{ns, N_DRTOn,s_nu,-1..m_mpc-1};              # Control variables mpc level
var ympc{ns, N_DRTOn,s_ny,1..p_mpc};                 # Output variables mpc level
var du_mpc{ns, N_DRTOn,s_nu,0..m_mpc-1};             # Delta_u = u(j,k+1) - u(j,k-1)            
var ysp_mpc{ns, N_DRTOn,s_ny,1..p_mpc};              # Set_point trajectory mpc      
var usp_mpc{ns, N_DRTOn,s_nu,0..m_mpc-1};            # Set_point trajectory 
var d_mpc{ns, N_DRTOn,s_ny};                         # Disturbance mpc


# Lagrange multiplier
var lambda1{ns, N_DRTOn,s_nx,0..p_mpc-1};             # Lambda equalities states
var lambda2{ns, N_DRTOn,s_ny,1..p_mpc};               # Lambda equalities outputs
var lambda3{ns, N_DRTOn,s_nu,0..m_mpc-1};             # Lambda equalities delta u           

#*********************************************************************************************************************
subject to

#---------------------------------------------------------------------------------------------------------------------
# DRTO Constraints
#---------------------------------------------------------------------------------------------------------------------
   
rto_x{s in ns, i in s_nx, j in 0..n_DRTO-1}:           xrto[s, i, j+1] = sum{n in s_nx} A_DRTO[s,i,n]*xrto[s, n, j] + 
                                                            sum{n in s_nu} B_DRTO[s, i, n]*urto[s, n, j];                 # Equivalent for x(k+1) = Ax(k)+Bu(k)

rto_y{s in ns, i in s_ny, j in 1..n_DRTO}:             yrto[s, i, j] = sum{n in s_nx} C_DRTO[s, i, n]*xrto[s, n, j] + d_drto[s,i];       # Equivalent for y(k) = Cx(k)+ d(0)

condition1{s in ns, i in s_nx}:                        xrto[s, i, 0] = xrto_0[s, i];                                            # Initialization

# Reference trajectory
# yref_hold1{s in ns, i in s_ny, j in 0..n_DRTO-sph by sph, jj in 2..sph}:    sph>=2 ==> yref[s, i, j+jj] = yref[s, i, j+1];    # Output Reference trajectory 1..n_DRTO
yref_hold2{s in ns, i in s_ny, j in n_DRTO + 1 .. n_DRTO + p_mpc -1}:         yref[s, i, j]    = yref[s, i, n_DRTO];          # Output Reference trajectory n_DRTO+1..n_DRTO + p_mpc-1

# uref_hold1{s in ns, i in s_nu, j in 0..n_DRTO-sph by sph, jj in 1..sph-1}: sph>=2 ==> uref[s, i, j+jj] = uref[s, i, j];     # Input Reference trajectory 0..n_DRTO-1
uref_hold2{s in ns, i in s_nu, j in n_DRTO..n_DRTO + m_mpc -2 }:              uref[s, i, j] = uref[s, i, n_DRTO-1];         # Input Reference trajectory n_DRTO..n_DRTO+m_mpc-2          


bounds_Y_rto{s in ns, i in s_ny, j in 1..n_DRTO}:              ymin[i] <= yrto[s, i, j] <= ymax[i];                                                        # Bounds on output variables

bounds_Y_ref{s in ns, i in s_ny, j in 1..n_DRTO+p_mpc-1}:      ymin_ref[i] <= yref[s, i, j]  <= ymax_ref[i];                                               # Bounds on output variables
bounds_Y_ref3{s in ns, i in s_ny}:                             sum {j in 1..n_DRTO} bref[s, i, j] <= N_ref;               # Bounds on output variables
bounds_Y_ref4{s in ns, i in s_ny, j in 1..n_DRTO-H_ref+1}:     sum {k in 0..H_ref-1} bref[s, i, j+k] <= 1;                         # Bounds on output variables
bounds_Y_ref5{s in ns, i in s_ny, j in 1..1}:                  bref[s, i, j] = 1;                                         # Bounds on output variables
bounds_Y_ref2{s in ns, i in s_ny, j in 1..n_DRTO-1}:           yref[s, i, j+1] - yref[s, i, j]  >= (ymin_ref[i] - ymax_ref[i])*bref[s, i, j+1];            # Bounds on output variables
bounds_Y_ref2a{s in ns, i in s_ny, j in 1..n_DRTO-1}:          yref[s, i, j+1] - yref[s, i, j]  <= (ymax_ref[i] - ymin_ref[i])*bref[s, i, j+1];            # Bounds on output variables

bound_U_ref{s in ns, i in s_nu,  j in 0..n_DRTO-1+m_mpc-1}:    umin_ref[i] <= uref[s, i, j] <= umax_ref[i];   # Bounds on control variables

##########################################################################################################################
# MPC SUBPROBLEMS FOR DRTO j = 0..n-1, where one MPC subproblem for every j
##########################################################################################################################

# ------------------------------------------------------------------------------------------------------------------------
# Linkage DRTP and MPC subproblems 
#-------------------------------------------------------------------------------------------------------------------------
# Disturbance model
condition_mpc_1{s in ns, n in 1..n_DRTOn-1, i in s_ny}:     d_mpc[s, n, i] = yrto[s, i, n] - sum{jj in s_nx} C_[i, jj]*xmpc[s, n-1, jj, 1];
condition_mpc_2{s in ns, i in s_ny}:                        d_mpc[s, 0, i] = ymeas[i]   - sum{jj in s_nx} C_[i, jj]*xmpc[s, 0, jj, 0];

# Control values u(-1) for each mpc subproblems
condition_mpc_3{s in ns, i in s_nu}:                        umpc[s, 0,i,-1] = urto_minus_one[i];      # previous inputs u for the FIRST KKT system, j = 0, k = 0
condition_mpc_4{s in ns, n in 1..n_DRTOn-1, i in s_nu}:     umpc[s, n,i,-1] = urto[s, i, n-1];        # previous inputs u for the REMAINING KKT system, j = 1..prto-1, k = 0

# Initialization states MPC subproblems
condition_mpc_6{s in ns, i in s_nx}:                        xmpc[s, 0,i,0]  = xmpc_0[i];               # x0 for first KKT subproblem N = 0
condition_mpc_7{s in ns, n in 1..n_DRTOn-1, i in s_nx}:     xmpc[s, n,i,0]  = xmpc[s, n-1, i, 1];      # x0 for remaining KKT subproblems n = 1..n_DRTO-1


# Set-points for mpc subproblems
YSPtraj_mpc1{s in ns, n in N_DRTOn, i in s_ny, j in 1..p_mpc}:          ysp_mpc[s, n, i, j] = yref[s, i, n+1];      # Set-point for output in mpc

USPtraj_mpc1{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:        usp_mpc[s, n, i, j] = uref[s, i, n];        # Set-point for inputs in mpc


#********************************************************************************************
# KKT contiditons for mpc subproblem
#*********************************************************************************************

#--------------------------------------------------------------------------------------------
# Primal Feasibility Equality constraints
#--------------------------------------------------------------------------------------------
# Linear state-space model constraints
primal1_mpc{s in ns, n in N_DRTOn, i in s_nx, j in 0..m_mpc-1}:
                    xmpc[s, n, i, j+1] = sum{jj in s_nx} A_[i, jj] * xmpc[s, n, jj, j] 
                                              +sum{jj in s_nu} B_[i, jj] * umpc[s, n, jj, j];

# Linear state-space model constraints
primal2_mpc{s in ns, n in N_DRTOn, i in s_nx, j in m_mpc..p_mpc-1}:
                    xmpc[s, n, i, j+1] = sum{jj in s_nx} A_[i, jj] * xmpc[s, n, jj, j] 
                                        +sum{jj in s_nu} B_[i, jj] * umpc[s, n, jj, m_mpc-1];

# Linear model for output constraints
primal3_mpc{s in ns, n in N_DRTOn, i in s_ny, j in 1..p_mpc}:
                    ympc[s, n, i, j] = sum {jj in s_nx} C_[i, jj] * xmpc[s, n, jj, j] 
                                                                 + d_mpc[s, n, i];

# Delta_u(k) = u(k)-u(k-1)
primal4_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:
                    du_mpc[s, n, i, j] = umpc[s, n, i, j] - umpc[s, n, i, j-1];


# Slack variables
var s_max{s in ns, N_DRTOn,s_nu,0..m_mpc-1} , >=0;         # Slack variables to convert inequality u<=u_max to equality
var s_min{s in ns, N_DRTOn,s_nu,0..m_mpc-1} , >=0;         # Slack variables to convert inequality u>=u_min to equality

# Multipliers for inequality
var mu_max{s in ns, N_DRTOn,s_nu,0..m_mpc-1} , >=0;        # Mu upper bound on u  s_min*mu1 = 0
var mu_min{s in ns, N_DRTOn,s_nu,0..m_mpc-1} , >=0;        # Mu lower bound on u  s_max*mu2 = 0 

# Complementarity conditions
var compl_cond = sum{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}  (mu_max[s, n, i, j] * s_max[s, n, i, j] + mu_min[s, n, i, j] * s_min[s, n, i, j]);

subject to

condition_mpc_5{s in ns, n in N_DRTOn, i in s_nu}:        urto[s, i, n]    = umpc[s, n, i, 0];            # feedback MPC-KKT u's to DRTO

# Inequalities expresses as equalities using slack variables u <= umax
primal5_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: umax[i] - umpc[s, n, i, j] - s_max[s, n, i, j] = 0; 
primal6_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: umpc[s, n, i, j] - umin[i] - s_min[s, n, i, j] = 0;

# Constraints on slack variables
primal7_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: s_max[s, n, i, j] >= 0;
primal8_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: s_min[s, n, i, j] >= 0;

#-----------------------------------------------------------------------------------------------------------------
# Dual feasibility 
#-----------------------------------------------------------------------------------------------------------------
dual1_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: mu_max[s, n, i, j] >= 0;
dual2_mpc{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}: mu_min[s, n, i, j] >= 0;

#--------------------------------------------------------------------------------------------------------------------
# Gradient of the Lagrangian
#----------------------------------------------------------------------------------------------------------------------

# Lagrange derivatives w.r.t y
lagrange_y{s in ns, n in N_DRTOn, i in s_ny, j in 1..p_mpc}: 2 * Q_[i] * ympc[s, n, i, j]
                                                            -2 * Q_[i] * ysp_mpc[s, n, i, j]
                                                            - lambda2[s, n,i,j] = 0;

# Lagrange derivatives w.r.t du
lagrange_du{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc - 1}: 2 * R_[i] * du_mpc[s, n, i, j] 
                                                        + lambda3[s, n, i, j] = 0;

# Lagrange derivatives w.r.t x
lagrange_xa{s in ns, n in N_DRTOn, i in s_nx, j in 0..p_mpc-2}: -lambda1[s, n, i, j]
                                                    + sum{jj in s_nx} A_[jj, i] * lambda1[s, n, jj, j+1]
                                                    + sum{jj in s_ny} C_[jj, i] * lambda2[s, n, jj, j+1] = 0;

lagrange_xb{s in ns, n in N_DRTOn, i in s_nx, j in p_mpc-1..p_mpc-1}:  -lambda1[s, n, i, j]
                                                            + sum{jj in s_ny} C_[jj, i] * lambda2[s, n, jj, j+1] = 0;
# Lagrange derivatives w.r.t u
lagrange_ua{s in ns, n in N_DRTOn, i in s_nu,  j in 0..m_mpc - 2}:  2 * S_[i] * umpc[s, n, i, j]
                                                        - 2 * S_[i] * usp_mpc[s, n, i, j]
                                                        + sum{jj in s_nx} B_[jj, i] * lambda1[s, n, jj, j]
                                                        - lambda3[s, n, i, j]
                                                        + lambda3[s, n, i, j+1] 
                                                        - mu_min[s, n, i, j]                                                       
                                                        + mu_max[s, n, i, j] = 0;

lagrange_ub{s in ns, n in N_DRTOn, i in s_nu, j in m_mpc-1..m_mpc-1}:       2 * S_[i] * umpc[s, n, i, j]
                                                                - 2 * S_[i] * usp_mpc[s, n, i, j]
                                                                + sum{k in m_mpc-1..p_mpc-1} sum{jj in s_nx} B_[jj, i] * lambda1[s, n, jj, k]
                                                                - lambda3[s, n, i, j]
                                                                - mu_min[s, n, i, j]
                                                                + mu_max[s, n, i, j] = 0;
                                                                
#--------------------------------------------------------------------------------------------------
# Bynary variables and BiG-M parameters to formulate complementariy constraints in non-linear form
#---------------------------------------------------------------------------------------------------
var cc_max{ns, N_DRTOn,s_nu,0..m_mpc-1} binary;        # For upper bound
var cc_min{ns, N_DRTOn,s_nu,0..m_mpc-1} binary;        # For lower bound
           
#_____________________________________________________________________________
# Constraints for linear and intenger formulation of complementariy conditons
#_____________________________________________________________________________
subject to

compl_max1{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:   s_max[s, n, i, j] - M_comp[i]*cc_max[s, n, i, j] <= 0;            # Upper bound complementarity conditions
compl_max2{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:   mu_max[s, n, i, j]- M_comp[i]*(1- cc_max[s, n, i, j]) <= 0;      # Upper bound complementarity conditions

compl_min1{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:   s_min[s, n, i, j] - M_comp[i]*cc_min[s, n, i, j] <= 0;        # Lower bound complementarity conditions
compl_min2{s in ns, n in N_DRTOn, i in s_nu, j in 0..m_mpc-1}:   mu_min[s, n, i, j] - M_comp[i]*(1-cc_min[s, n, i, j]) <= 0;   # Lower bound complementarity conditions

##################################################################################################################
# Scheduling variables and parameters
################################################################################################################## 
var   ysp{ns, i in s_ny,0..n_DRTO};           # Target set point for DRTO
param inventory_0{grades};                    # Quantity of each product grade in inventory at time zero
param molar_mass{grades};                     # Molar mass
param ytargets{grades, s_ny};                 # Output targets for every grade
param epsilon{grades,  s_ny};                 # For tolerance band on set-point

#__________________________________________________________________________________________________________
# Defining the grade being produced
#__________________________________________________________________________________________________________
var binary_grade{ns, 1..3, grades, i in s_ny, 1..n_DRTO} binary;
var gamma_tilda{ns, grades,0..n_DRTO} binary;
var beta_tilda{ns, grades,0..n_DRTO} >=0, <=1;
var mass_g{ns, grades, 1..n_DRTO};
var gammai{ns, grades, 0..n_DRTO} binary;
param H_z integer > 0;
param N_z integer > 0;
param Lz  integer > 0;
param zplant{grades, -Lz+1..0} binary;
var z_targs{ns, grades,-Lz+1..n_DRTO} >=0, <=1;

subject to
plant_info1{s in ns, g in grades, j in -Lz+1..0}:                              z_targs[s, g, j] - zplant[g, j] <=  (1 - gamma_tilda[s, g, 1 + H_z-1]);
plant_info2{s in ns, g in grades, j in -Lz+1..0}:                              z_targs[s, g, j] - zplant[g, j] >= -(1 - gamma_tilda[s, g, 1 + H_z-1]);
plant_info3{s in ns, g in grades, j in -Lz+1..0}:                              z_targs[s, g, j] <= zplant[g, j];

plant_info4{s in ns, g in grades}:   sum{k in -H_z+1..0} zplant[g, k] = H_z  ==>  gammai[s, g, 0] = 1;

one_product_per_step1{s in  ns, n in 0..n_DRTO}:                            sum{g in grades} gammai[s, g, n] = 1;
define_target{s in ns, i in s_ny, n in 0..n_DRTO}:                          ysp[s, i,n] = sum{g in grades} gammai[s, g, n]*ytargets[g,i];

one_product_per_step2{s in ns, n in 1..n_DRTO, i in s_ny, g in grades}:     sum{l in 1..3} binary_grade[s, l, g, i, n] = 1;  

within_target1{s in ns, n in 1..n_DRTO, g in grades}:                       sum{i in s_ny} binary_grade[s, 2, g, i, n] <= ny_ - (1 - z_targs[s, g, n]);    
within_target2{s in ns, n in 1..n_DRTO, g in grades}:                       sum{i in s_ny} binary_grade[s, 2, g, i, n] >= ny_*z_targs[s, g, n];  

#exclude_poss1{s in ns, n in 1..n_DRTO, g in grades}:                        z_targs[s, g, n] <= gammai[s, g, n];
exclude_poss2{s in ns, n in 1..n_DRTO}:                                     sum{g in grades} z_targs[s, g, n] <= 1;

z_hold1{s in ns, g in grades, n in - H_z + 1..n_DRTO - H_z + 1}:            sum{k in 0..H_z-1} z_targs[s,g,n+k] <= H_z - (1 - gamma_tilda[s, g, n + H_z-1]);
z_hold2{s in ns, g in grades, n in - H_z + 1..n_DRTO - H_z + 1}:            sum{k in 0..H_z-1} z_targs[s,g,n+k] >= H_z*gamma_tilda[s, g, n + H_z-1];
z_hold3{s in ns, g in grades, n in 0..n_DRTO}:                              sum{k in 0..Lz-1}  z_targs[s,g,n-k] >= Lz*beta_tilda[s, g, n];
z_hold4{s in ns, g in grades, n in 1..n_DRTO}:                              gamma_tilda[s,g,n] <= z_targs[s,g,n];

store_1{s in ns, g in grades}:                                              sum{n in 0..n_DRTO} beta_tilda[s, g, n]        <= N_z - gammai[s, g, n_DRTO];                             # Can only produce a given product once in the scheduling horizon                      
store_2{s in ns, n in 0..n_DRTO-1, g in grades}:                            gammai[s, g, n+1]  - gammai[s, g, n]           >= -beta_tilda[s, g, n];
store_3{s in ns, n in 0..n_DRTO-1, g in grades}:                            gamma_tilda[s, g, n+1] - gamma_tilda[s, g, n]  >= -beta_tilda[s, g, n];
store_4{s in ns, g in grades, n in 0..n_DRTO-1}:                            gammai[s, g, n+1]                              <=  1 - beta_tilda[s, g, n];
store_5{s in ns, g in grades, n in 1..n_DRTO}:                              beta_tilda[s, g, n]                            <=  gamma_tilda[s, g, n];
store_6{s in ns, g in grades, n in 1..n_DRTO}:                              gammai[s, g, n]                                >=  gamma_tilda[s, g, n];
store_7{s in ns, n in 0..n_DRTO}:                                           sum{g in grades} beta_tilda[s, g, n]           <= 1;

#_________________________________________________
# Target band
#_________________________________________________
subject to

lower_target_band_yrto{s in ns, i in s_ny, g in grades, n in 1..n_DRTO}:   yrto[s, i, n] >= binary_grade[s, 1, g, i, n] * ymin[i] 
                                                                                          + binary_grade[s, 2, g, i, n] * (ytargets[g, i] - epsilon[g, i])
                                                                                          + binary_grade[s, 3, g, i, n] * (ytargets[g, i] + epsilon[g, i]);  

upper_target_band_yrto{s in ns, i in s_ny, g in grades, n in 1..n_DRTO}:   yrto[s, i, n] <= binary_grade[s, 1, g, i, n] * (ytargets[g, i] - epsilon[g, i]) 
                                                                                          + binary_grade[s, 2, g, i, n] * (ytargets[g, i] + epsilon[g, i])
                                                                                          + binary_grade[s, 3, g, i, n] * ymax[i];


mass_constraint1{s in ns, g in grades, n in 1..n_DRTO}:    mass_g[s, g, n] >= -gamma_tilda[s, g, n] * 120;
mass_constraint2{s in ns, g in grades, n in 1..n_DRTO}:    mass_g[s, g, n] <= gamma_tilda[s, g, n] * 120;
mass_constraint3{s in ns, g in grades, n in 1..n_DRTO}:    mass_g[s, g, n] - (yrto[s, 1, n] + lin_point_yrto[1])*flow_rate*molar_mass[g]*dt_MPC >=   - (1 - gamma_tilda[s, g, n]) * 120;
mass_constraint4{s in ns, g in grades, n in 1..n_DRTO}:    mass_g[s, g, n] - (yrto[s, 1, n] + lin_point_yrto[1])*flow_rate*molar_mass[g]*dt_MPC <=     (1 - gamma_tilda[s, g, n]) * 120;


#______________________________________________________
# Inventory
#______________________________________________________ 
param demand_nominal{ns, grades, 1..n_DRTO} default 0;            # Quantity demanded of each grade at every time step
param lb_inv{grades};                                             # Lower bound on inventory
param ub_inv{grades};                                             # Upper bound on inventory
var inventory{ns, grades, 0..n_DRTO}, >=0;                        # Amount of material in inventory
var demand{s in ns, g in grades, n in 1..n_DRTO}, >=0, <= demand_nominal[s, g, n];                  # Demand Met
var demand_not_met{s in ns, g in grades, n in 1..n_DRTO}, >=0, <= demand_nominal[s, g, n];          # Fraction of the demand that is met; set to zero to enforce that all the demand is met

subject to

demand_con2{s in ns, g in grades, n in 1..n_DRTO}:               demand[s, g, n] = demand_nominal[s, g, n] - demand_not_met[s, g, n];
inventory_constraint1{s in ns, g in grades}:                     inventory[s, g, 0] = inventory_0[g];                                                                 # Initial amount of product in inventory
inventory_constraint2{s in ns, g in grades, n in 1..n_DRTO}:     inventory[s, g, n] = inventory[s, g, n-1]     
                                                                                        + mass_g[s, g, n]     
                                                                                        - demand[s, g, n];

#--------------------------------------------------------------------------------------------
# OBJECTIVE FUNCTION
#--------------------------------------------------------------------------------------------
param cost_inventory{grades};                 # Storage cost for each product grade
param cost_backorder{grades};                 # Cost of not meeting the demand
param cost_input{1..nu_};                     # Cost of the input
param cost_product{grades};                   # Cost final product
param value_inv{grades};                      # Value inventory

var cost_raw_material{s in ns} =  1/100*(sum{n in 0..n_DRTO-1} cost_input[1]*(urto[s, 1,n]+lin_point_urto[1])*flow_rate*dt_MPC);
var cost_inv_g {ns, grades}; 
var value_inv_g {ns, grades} >= 0;
var cost_backorder_g {ns, grades}; 
var revenue {ns, grades};
var profit_sce{ns};


maximize

economic_objective:  sum{s in ns} prob_scenario[s]*(-cost_raw_material[s] + profit_sce[s]);


subject to

con_obj1{s in ns, g in grades}: cost_backorder_g[s, g] >= sum{n in 1..n_DRTO} cost_backorder[g]*demand_not_met[s, g, n]; 
con_obj2{s in ns, g in grades}: revenue[s, g] <= sum{n in 1..n_DRTO} cost_product[g]*demand[s, g, n];
con_obj3{s in ns, g in grades}: cost_inv_g[s, g] >= sum{n in 1..n_DRTO} cost_inventory[g]*inventory[s, g, n-1];
con_obj4{s in ns}:              profit_sce[s] <= 1/100*(sum{g in grades} (revenue[s, g] + value_inv_g[s, g]- cost_inv_g[s, g] - cost_backorder_g[s, g]));
con_obj5{s in ns, g in grades}: value_inv_g[s, g] <= value_inv[g]*inventory[s, g, n_DRTO];
con_obj6{s in ns, g in grades}: value_inv_g[s, g] <= value_inv[g]*ub_inv[g];



#----------------------------------------------------------------------------------------------------------------------------
# NON-ANTICIPATIVITY CONSTRAINTS
#----------------------------------------------------------------------------------------------------------------------------
param n_period integer > 0;       # Number of periods
param tD{0..n_period};            # Time-step of uncertainty realization   
param n_bundles{1..n_period};                                 # Number of bundles in each period
set  bundles{i in 1..n_period, 1..n_bundles[i]} ordered;      # Bundles

subject to

NAC1yref{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], i in s_ny, n in tD[p-1]+1..tD[p]}:         yref[member(1, bundles[p, j]), i, n] = yref[s, i, n];

NAC1urto{p in 1..1, j in 1..n_bundles[p], s in bundles[p, j], i in s_nu, n in 1..1}:                            urto[member(1, bundles[p, j]), i, n-1] = urto[s, i, n-1];

#NAC1xrto{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], i in s_nx, n in tD[p-1]+1..tD[p]}:        xrto[member(1, bundles[p, j]), i, n] = xrto[s, i, n];

#NAC1m{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], g in grades, n in tD[p-1]+1..tD[p]}:         mass_g[member(1, bundles[p, j]), g, n] = mass_g[s, g, n];

NAC1bref{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], i in s_nu, n in  tD[p-1]+1..tD[p]}:        bref[member(1, bundles[p, j]), i, n]   = bref[s, i, n];

#NAC1bg{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], k in 1..3, i in s_ny, g in grades, n in tD[p-1]+1..tD[p]}:         binary_grade[member(1, bundles[p, j]), k, g, i, n] = binary_grade[s, k, g, i, n];

#NAC2bt{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], g in grades, n in tD[p-1]+1..tD[p]}:        binary_tilda[member(1, bundles[p, j]), g, n] = binary_tilda[s, g, n];

#NAC2gamma{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], g in grades, n in tD[p-1]+1..tD[p]}:     gammai[member(1, bundles[p, j]), g, n] = gammai[s, g, n];

#NAC2ztarg{p in 1..n_period, j in 1..n_bundles[p], s in bundles[p, j], g in grades, n in tD[p-1]+1..tD[p]}:     z_targs[member(1, bundles[p, j]), g, n] = z_targs[s, g, n];

NAC2cc_max{p in 1..1, j in 1..n_bundles[p], s in bundles[p, j], i in s_nu, n in 1..1, k in 0..m_mpc-1}:         cc_max[member(1, bundles[p, j]), n-1, i, k] = cc_max[s, n-1, i, k];

NAC2cc_min{p in 1..1, j in 1..n_bundles[p], s in bundles[p, j], i in s_nu, n in 1..1, k in 0..m_mpc-1}:         cc_min[member(1,bundles[p, j]), n-1, i, k] = cc_min[s, n-1, i, k];


#lower_target_band_ref{s in ns, i in s_ny, n in tD[1]+1..n_DRTO}:     yref[s, i, n] >= (sum{g in grades} gamma_tilda[s, g, n]*(ytargets[g, i] - 0*epsilon[g, i])) + (1 -  sum{g in grades} gamma_tilda[s, g, n]) * ymin_ref[i];

#upper_target_band_ref{s in ns, i in s_ny, n in tD[1]+1..n_DRTO}:     yref[s, i, n] <= (sum{g in grades} gamma_tilda[s, g, n]*(ytargets[g, i] + 0*epsilon[g, i])) + (1 -  sum{g in grades} gamma_tilda[s, g, n]) * ymax_ref[i];


param total_solve_time = _total_solve_elapsed_time;   # elapsed time for all solve commands
param solve_time = _solve_elapsed_time;               # elapsed seconds for most recent solve command;    
#param compl_cond default 0;                          # To check complementarity conditons
param obj_sce{s in ns} default 0;                     # Objective values for each scenario
option presolve_eps 1e-7;

    