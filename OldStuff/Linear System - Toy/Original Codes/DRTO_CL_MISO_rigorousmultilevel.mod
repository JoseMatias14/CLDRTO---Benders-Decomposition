#===============================================================================
# Closed-loop DRTO (linear) model with embedded linear MPC (KKT constraints)
# Mohammad Zamry Jamaludin <jamalumz@mcmaster.ca>
# Date: 03/04/2014
# Model: MISO transfer function
#-------------------------------------------------------------------------------

# User-defined parameters
#-------------------------------------------------------------------------------
# RTO parameters
param prto_;               # prediction horizon of RTO 

# Model parameters
param nx_;                 # No. of differential states
param ny_;                 # No. of outputs
param nu_;                 # No. of inputs
param A_{1..nx_, 1..nx_};  # A matrix
param B_{1..nx_, 1..nu_};  # B matrix
param C_{1..ny_, 1..nx_};  # C matrix

# Tuning parameters
param Q_{1..ny_};          # Q matrix (diagonal elements only)
param R_{1..nu_};          # R matrix (diagonal elements only)
param S_{1..nu_};          # S matrix (diagonal elements only)
param p_;                  # prediction horizon
param m_;                  # control horizon
param sph_;                # set-point hold

# Time parameters and variable sets
set J_  := 0..prto_ - 1;   # RTO prediction
set P_  := 0..p_ - 1;      # prediction
set M_  := 0..m_ - 1;      # control
set I_  := 1..2;           # inequality hard constraints wrt umin & umax
set Ix_ := 1..nx_;         # state variables
set Iy_ := 1..ny_;         # output variables
set Iu_ := 1..nu_;         # input variables

# Optimization bounds
param xl_{Ix_};            # lower bound of state var in RTO
param xu_{Ix_};            # upper bound of state var in RTO
param yl_{Iy_};            # lower bound of output var in RTO
param yu_{Iy_};            # upper bound of output var in RTO
param umin_{1..nu_};       # minimum u in MPC
param umax_{1..nu_};       # maximum u in MPC
param dumin_{Iu_};         # lower bound of du in MPC
param dumax_{Iu_};         # upper bound of du in MPC
param xmin_{Ix_};          # lower bound of state var in MPC
param xmax_{Ix_};          # upper bound of state var in MPC
param ymin_{Iy_};          # lower bound of output var in MPC
param ymax_{Iy_};          # upper bound of output var in MPC
param yspmin_{Iy_};        # lower bound of ysp's in MPC
param yspmax_{Iy_};        # upper bound of ysp's in MPC
param uspmin_{Iu_};        # lower bound of usp's in MPC
param uspmax_{Iu_};        # upper bound of usp's in MPC
param d_{Iy_, 1..prto_};   # constant output disturbance estimate

# Computed sets, parameters and variables
#-------------------------------------------------------------------------------

# Initialization of DRTO execution
param xrtoinit_{Ix_};		       # initial condition of x
param yrtotarget_{Iy_, 1..prto_};  # output target
param urtoprev_{Iu_};        	   # previous plant input


# Linking
param xss_{Ix_};                   #  steady-state values of state variables
param yss_{Iy_};                   #  steady-state values of output variables
param uss_{Iu_};                   #  steady-state values of input variables

# KKT variables
var x_   {xkkt in Ix_, J_, 0..p_}	   ;#, >=  xmin_[xkkt]  , <= xmax_[xkkt];      # x_kkt vector
var y_   {ykkt in Iy_, J_, 1..p_}      ;#, >=  ymin_[ykkt]  , <= ymax_[ykkt];      # y_kkt vector
var yref_{yref  in Iy_, J_, 1..p_}     ;#, >=  yspmin_[yref] , <= yspmax_[yref];   # ysp in KKT system
var uref_{uref  in Iu_, J_, 0..m_ - 1} ;#, >=  uspmin_[uref] , <= uspmax_[uref];   # usp in KKT system
var du_{du in Iu_, J_, 0..m_ - 1}      ;#, >=  dumin_[du] , <= dumax_[du];	       # du_kkt vector
var u_{Iu_, J_, -1..m_ - 1};	                                                   # u_kkt vector 
var s1_{Iu_, J_, 0..m_ - 1}            , >=  0;    	                               # upper slack inequalities
var s2_{Iu_, J_, 0..m_ - 1}            , >=  0;    	                           	   # lower slack inequalities

var lamb1_{Ix_, J_, 1..p_}             ;#, >= -1e9 , <= 1e9;   		               # lagrange multipliers for state predictions
var lamb2_{Iy_, J_, 1..p_}             ;#, >= -1e9 , <= 1e9;   		               # lagrange multipliers for output predictions
var lamb3_{Iu_, J_, 0..m_ - 1}         ;#, >= -1e9 , <= 1e9;   		               # lagrange multipliers for control inputs 0..p-1
var mu1_{Iu_, J_, 0..m_ - 1}            , >=  0;    	                               # lagrange multipliers for upper slack inequalities
var mu2_{Iu_, J_, 0..m_ - 1}            , >=  0;    	                               # lagrange multipliers for lower slack inequalities

# DRTO variables
var ysp_{ysp in Iy_, 0..prto_ - 1 + p_ - 1}    , >=  yspmin_[ysp] , <= yspmax_[ysp];     # DRTO set-point y for each KKT system
var usp_{usp in Iu_, 0..prto_ - 1 + m_ - 1}    , >=  uspmin_[usp] , <= uspmax_[usp];     # DRTO set-point u for each KKT system
var xrto_{xrto in Ix_, 0..prto_}               , >=  xl_[xrto] , <= xu_[xrto];	         # x_rto vector
var yrto_{yrto in Iy_, 1..prto_}               , >=  yl_[yrto] , <= yu_[yrto];	         # y_rto vector (Case: SP single step change)
var urto_{Iu_, 0..prto_ - 1};							                                 # u_rto vector (input constraints imposed on the inner subproblems - KKT)
 

subject to

#-------------------------------------------------------------------------------
# RTO closed-loop prediction
# Remark: A, B, C state-space matrices are written in deviation form.
# At the RTO level, the model is not in deviation form (actual variable values.
# Deviation-form model can be obtained appropriately through substraction from 
# the steady-states, e.g.(xrto - xss).
#-------------------------------------------------------------------------------
rto0_{i_ in Ix_, j_ in 0..prto_ - 1}:                    (xrto_[i_, j_ + 1] - xss_[i_]) =   sum{f_ in Ix_}(A_[i_, f_] * (xrto_[f_, j_] - xss_[f_]))
                                                                                          + sum{f_ in Iu_}(B_[i_, f_] * (urto_[f_, j_] - uss_[f_]));
rto1_{i_ in Iy_, j_ in 1..prto_}:                        (yrto_[i_, j_] - yss_[i_])     =   sum{f_ in Ix_}(C_[i_, f_] * (xrto_[f_, j_] - xss_[f_]));
#rto2_{i_ in Iy_, j_ in 1..prto_}:                        (yrto_[i_, j_] - yss_[i_])     =   sum{f_ in Ix_}(C_[i_, f_] * (xrto_[f_, j_] - xss_[f_])) + d_[i_, j_];

#-------------------------------------------------------------------------------
# Linkages between RTO and MPC
#-------------------------------------------------------------------------------
# Initialization for DTO
condition1_{i_ in Ix_}:                                  xrto_[i_, 0]   = xrtoinit_[i_];                # initial states x for RTO system, j = 0

# Initialization for KKT, and linkages between DRTO & KKT
condition2_{i_ in Ix_}:                                  x_[i_, 0, 0]   = xrtoinit_[i_] - xss_[i_];     # initial states x for the FIRST KKT system, j = 0, k = 0
condition3_{i_ in Ix_, j_ in 1..prto_ - 1}:              x_[i_, j_, 0]  = xrto_[i_, j_] - xss_[i_];     # initial states for the REMAINING KKT system, j = 1..prto-1, k = 0
condition4_{i_ in Iu_}:                                  u_[i_, 0, -1]  = urtoprev_[i_] - uss_[i_];     # previous inputs u for the FIRST KKT system, j = 0, k = 0
condition5_{i_ in Iu_, j_ in 1..prto_ - 1}:              u_[i_, j_, -1] = urto_[i_, j_ - 1] - uss_[i_]; # previous inputs u for the REMAINING KKT system, j = 1..prto-1, k = 0
condition6_{i_ in Iu_, j_ in 0..prto_ - 1}:              urto_[i_, j_]  = u_[i_, j_, 0] + uss_[i_];     # feedback MPC-KKT u's to RTO

#-------------------------------------------------------------------------------
# Relationship between yref and ysp:= shifted reference trajectory
#-------------------------------------------------------------------------------
# CV set-points
YSPtraj1_{i_ in Iy_, j_ in J_, k_ in 1..p_}:                                 yref_[i_, j_, k_] = ysp_[i_, j_ + k_ - 1];
YSPtraj2_{i_ in Iy_, j_ in prto_ + 1..prto_ - 1 + p_ - 1}:                   ysp_[i_, j_] = ysp_[i_, j_ - 1];
YSPhold1_{i_ in Iy_, j_ in prto_ - sph_..prto_ - 1 + p_ - 1}:                ysp_[i_, j_] = ysp_[i_, j_ - 1];
YSPhold2_{i_ in Iy_, j_ in 0..prto_ - 2 by sph_, jj_ in 0..sph_ - 2}:	     ysp_[i_, j_ + jj_] = ysp_[i_, j_ + jj_ + 1]; 

# To comment a section: /* ..... */
# MV set-points
USPtraj1_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:              				 uref_[i_, j_, k_] = usp_[i_, j_ + k_];
USPtraj2_{i_ in Iu_, j_ in prto_ + 1..prto_ - 1 + m_ - 1}:    				 usp_[i_, j_] = usp_[i_, j_ - 1];
USPhold1_{i_ in Iu_, j_ in prto_ - sph_..prto_ - 1 + m_ - 1}:                usp_[i_, j_] = usp_[i_, j_ - 1];
USPhold2_{i_ in Iu_, j_ in 0..prto_ - 2 by sph_, jj_ in 0..sph_ - 2}:	     usp_[i_, j_ + jj_] = usp_[i_, j_ + jj_ + 1]; 

#-------------------------------------------------------------------------------
# End-point constraints
#-------------------------------------------------------------------------------
yend1_{i_ in Iy_}:		                                   ysp_[i_, prto_ - 1] = yrtotarget_[i_, prto_];
yend2_{i_ in Iy_}:                                         ysp_[i_, prto_ - 1] = yrto_[i_, prto_] ;
#yend3_{i_ in Iy_, j_ in prto_ - sph_..prto_}:             yrto_[i_, j_] = yrto_[i_, j_ - 1];
uend1_:		                                               usp_[2, prto_ - 1]  = urto_[2, prto_ - 1];
#uend2_{i_ in Iu_, j_ in prto_ - sph_..prto_ - 1}:		   urto_[i_, j_] = urto_[i_, j_ - 1];
#uend3_{j_ in prto_ - sph_..prto_ - 1}:		               urto_[2, j_] = urto_[2, j_ - 1] ;

#-------------------------------------------------------------------------------
# KKT systems for MPC's QP subproblems
#-------------------------------------------------------------------------------
# Primal feasibility - equality constraints
primal1_{i_ in Ix_, j_ in J_, k_ in 0..m_ - 1}:            x_[i_, j_, k_ + 1] =   sum{f_ in Ix_}(A_[i_, f_] * x_[f_, j_, k_]) 
                                                                                + sum{f_ in Iu_}(B_[i_, f_] * u_[f_, j_, k_]);        		# state prediction

primal2_{i_ in Ix_, j_ in J_, k_ in m_..p_ - 1}:           x_[i_, j_, k_ + 1] =   sum{f_ in Ix_}(A_[i_, f_] * x_[f_, j_, k_]) 
                                                                                + sum{f_ in Iu_}(B_[i_, f_] * u_[f_, j_, m_ - 1 ]);  		# state prediction
																			  
primal3_{i_ in Iy_, j_ in J_, k_ in 1..p_}:                y_[i_, j_, k_] = sum{f_ in Ix_}(C_[i_, f_] * x_[f_, j_, k_]) + d_[i_, j_ + 1];   # output prediction

primal4_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:            du_[i_, j_, k_] = u_[i_, j_, k_] - u_[i_, j_, k_ - 1];                     		# input change 0..p-1

# Primal feasibility - inequality constraints expressed as equalities using slack variables
primal5_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:            u_[i_, j_, k_] - umin_[i_] - s1_[i_, j_, k_] = 0;   # umin bounds
primal6_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:            umax_[i_] - u_[i_, j_, k_] - s2_[i_, j_, k_] = 0;   # umax bounds

# Primal feasibility - inequality hard constraints
primal7_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:		       s1_[i_, j_, k_] >= 0;
primal8_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:		       s2_[i_, j_, k_] >= 0;

# Dual feasibility
dual1_  {i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:		       mu1_[i_, j_, k_] >= 0; 
dual2_  {i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:		       mu2_[i_, j_, k_] >= 0; 

# Lagrange derivative w.r.t y
lagrange1_{i_ in Iy_, j_ in J_, k_ in 1..p_}:              2 *    y_[i_, j_, k_] * Q_[i_] 
                                                         - 2 * yref_[i_, j_, k_] * Q_[i_] 
                                                         - lamb2_[i_, j_, k_] = 0;
														 
# Lagrange derivative w.r.t du
lagrange2_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1}:          2 * du_[i_, j_, k_] * R_[i_] 
                                                         + lamb3_[i_, j_, k_] = 0;
														 
# Lagrange derivative w.r.t x
lagrange3_{i_ in Ix_, j_ in J_, k_ in 0..p_ - 2}:        - lamb1_[i_, j_, k_ + 1]
                                                         + sum{f_ in Ix_}(A_[f_, i_] * lamb1_[f_, j_, k_ + 2])
                                                         + sum{f_ in Iy_}(C_[f_, i_] * lamb2_[f_, j_, k_ + 1]) = 0;
													 											 
lagrange4_{i_ in Ix_, j_ in J_, k_ in p_ - 1..p_ - 1}:   - lamb1_[i_, j_, k_ + 1] 
                                                         + sum{f_ in Iy_}(C_[f_, i_] * lamb2_[f_, j_, k_ + 1]) = 0;
														 
# Lagrange derivative w.r.t u
lagrange5_{i_ in Iu_, j_ in J_, k_ in 0..m_ - 2}:          2 *    u_[i_, j_, k_] * S_[i_]
														 - 2 * uref_[i_, j_, k_] * S_[i_]
														 + sum{f_ in Ix_}(B_[f_, i_] * lamb1_[f_, j_, k_ + 1])
                                                         - lamb3_[i_, j_, k_]
                                                         + lamb3_[i_, j_, k_ + 1]
                                                         - mu1_[i_, j_, k_]
                                                         + mu2_[i_, j_, k_] = 0;
														 
lagrang6_{i_ in Iu_, j_ in J_, k_ in m_ - 1.. m_ - 1}:     2 *    u_[i_, j_, k_] * S_[i_]
														 - 2 * uref_[i_, j_, k_] * S_[i_]    
														 + sum{kk_ in m_ - 1..p_ - 1} sum{f_ in Ix_}(B_[f_, i_] * lamb1_[f_, j_, kk_ + 1]) 
                                                         - lamb3_[i_, j_, k_]
                                                         - mu1_[i_, j_, k_]
                                                         + mu2_[i_, j_, k_] = 0;
												 

#------------------------------------------------------------------------------
# Load data from MATLAB
#------------------------------------------------------------------------------
data;
include matlabout.dat;
include initialization.dat;
#------------------------------------------------------------------------------

# === Objective function === 
minimize obj:
#------------------------------------------------------------------------------
# Economics or target tracking
#------------------------------------------------------------------------------ 
#1e1 * sum{i_ in Iy_, j_ in 1..prto_}  ( yrto_[i_, j_] - yrtotarget_[i_, j_]   )^2 +

- 1e2* sum{j_ in 1..prto_} ( ( 1 ) *
                           ( 0.5*tanh( 20*(yrto_[1, j_] - yrtotarget_[1, j_] + 0.1*yrtotarget_[1, j_]) ) + 0.5 )^1 *   # lb econ
                           ( 0.5*tanh( 20*(yrtotarget_[1, j_] + 0.1*yrtotarget_[1, j_] - yrto_[1, j_]) ) + 0.5 )^1 )   # up econ

+ 1e0 * sum{j_ in 0..prto_ - 1}  ( 2*urto_[1, j_] + 10*urto_[2, j_]  )  						   
#------------------------------------------------------------------------------
# Complementarity constraints - exact penalty formulation
+ 80 * sum{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1} ( mu1_[i_, j_, k_] * s1_[i_, j_, k_] + 
                                                  mu2_[i_, j_, k_] * s2_[i_, j_, k_] );
											   
#------------------------------------------------------------------------------
# Solver settings
#------------------------------------------------------------------------------
# --- Selected solver ---
#solexpand > RTOsiso_bj.txt;
option ipopt_options "mu_init=1e-1   acceptable_tol=1e-9   max_iter=3000   linear_solver=ma27";
#option ipopt_options "mu_init=1e-2   acceptable_tol=1e-9";
option solver ipopt;
#option solver conopt;
# -----------------------

# --- AMPL-specific ---
#option presolve 0;
option eexit 0;
option presolve_fixeps 0;
option presolve_warnings 100;
option show_stats 1;
# ---------------------

# --- AMPL optional settings ---
# Turn on if: you want to substitute all variables and reduce the system
# option substout 1;
# --- CONOPT-specific ---
# --- KNITRO-specific ---
# option knitro_options 'outlev=2';
# --- IPOPTC-specific ---
# if match($solver, 'ipoptc') != 0 then { option ipopt_options 'iprint=0'; }
#--------------------------------

# --- Solve the problem ---
solve;
#display x_, y_, u_, s_, mu_, lamb0_, lamb1_, lamb2_, lamb3_, lamb4_, lamb6_, lamb7_, xbar_, ubar_, ysp_, xrto_, yrto_, urto_ >> PSinitguess.txt; # display and save variables 
#---------------------------


#------------------------------------------------------------------------------
# Check Cost
#------------------------------------------------------------------------------
param Cost_;
let Cost_ := sum{j_ in 0..prto_ - 1}  ( 0*urto_[1, j_] + 1*urto_[2, j_]  )    ;										
printf "Input Cost: %.4f\n", Cost_;

#------------------------------------------------------------------------------
# Check final complementarity value
#------------------------------------------------------------------------------
param complementvalue_;
let complementvalue_ := sum{i_ in Iu_, j_ in J_, k_ in 0..m_ - 1} ( mu1_[i_, j_, k_] * s1_[i_, j_, k_] + 
                                                                    mu2_[i_, j_, k_] * s2_[i_, j_, k_] );										
printf "Final complementarity value: %.15f\n", complementvalue_;
check complementvalue_ <= 1e-09;
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Output to MATLAB
#------------------------------------------------------------------------------
include amplconduits.inc;
#shell "amplconduits.py";
#------------------------------------------------------------------------------

display complementvalue_;



