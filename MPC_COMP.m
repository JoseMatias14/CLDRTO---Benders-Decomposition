function [ur_k,yr_k,V_k]=MPC_COMP(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ad,Bd,Cd,Ap,Bp,Cp,umax,umin,dumax,x0,u0,ys,distModel)
%  Simulates the closed-loop system with MPC with state-space model i
%  ur,yr,Vk  - Input and output responses (dimension: nu x nsim | ny x
%  nsim) and OF trajectory
%  p    - Optimization horizon
%  m    - Control horizon
%  nu   - Number of inputs
%  ny   - Number of outputs
%  nx   - Dimension of the state vector
%  nsim - Simulation time
%  q  - Output weights (dimension: 1 x ny)
%  r -  Input weights (dimension: 1 x nu)
%  A,B,C - State, input and output matrices of the state-space model used in the MPC controller
%  Ad,Bd,Cd - State, input and output matrices of the state-space model
%  used in the MPC controller (with disturbance model)
%  Ap,Bp,Cp - State, input and output matrices of the state-space model used to represent the true plant
%  umax,umin - Max and min values for the inputs (dimension: ny x 1)
%  dumax - Max input change (dimension: ny x 1)
%  x0 - initial state value
%  u0 - initial input value
%  ys   - Set-points for the outputs (dimension: ny x 1)
% distModel - flag to indicate if dist model is used (1) or not (0)
%               if flag is 2, we compute an unconstrained MPC

%% cheking if disturbance model is used
if distModel ~= 0
    A_mpc = Ad;
    B_mpc = Bd;
    C_mpc = Cd;
    nx_mpc = nx + ny;
else 
    A_mpc = A;
    B_mpc = B;
    C_mpc = C;
    nx_mpc = nx;
end

%% Building MPC matrices (PQI 5780 – Chemical Processes Control I - pg. 16-18)
Phi=[]; % matrix multiplying  x(k)
ThA=[]; % temp matrix for creating matrix that multiplies u array
for in=1:p
    Phi=[Phi;C_mpc*A_mpc^in];
    ThA=[ThA;C_mpc*A_mpc^(in-1)*B_mpc];
end

% Creating the Dynamic Matrix (that multiplies u array)
a = ThA;
DynM =a;

%N.B. code is only valid if m >= 2
for iu = 1:m - 2
    a = [zeros(ny,nu);a(1:(p-1)*ny,:)];
    DynM=[DynM a];
end

% adjusting dynamic matrix for since p > m (last column)
b = C_mpc*B_mpc;
Ai = eye(nx_mpc);
for in = 1:p-m
    Ai = Ai + A_mpc^in;
    b = [b;C_mpc*Ai*B_mpc];
end
Theta=[DynM [zeros(ny*(m-1),nu);b]];

% Creating Qbar and Rbar matrices
Qbar = diag(repmat(q,[1,p]));
Rbar = diag(repmat(r,[1,m]));

% Creating input movement OF penalty matrix 
M=[zeros((m-1)*nu,nu) eye(nu*(m-1));zeros(nu) zeros(nu,nu*(m-1))];
Ibar=[eye(nu);zeros(nu*(m-1),nu)];
IM = eye(nu*m) - M';

%Matrix H
H = Theta'*Qbar*Theta+IM'*Rbar*IM;

% Creating setpoint array (ny X nsim)
ysp= repmat(ys,[p,1]);

% Creating constraints array (ny X nsim)
Dumax = repmat(dumax,[m,1]);
Umax = repmat(umax,[m,1]);
Umin = repmat(umin,[m,1]);

%% State observer
Kf = FKalman(ny,A,C,100);

%% Starting simulation
% model and plant state (initial state is known)
xmk = x0;
xpk = xmk;

% model and plant outputs 
ymk = C*xmk;
ypk = Cp*xpk;
% prediction error
de = ypk - ymk;

% initial input (for computing input movement du)
uk_1 = u0;

V_k = [];
ur_k = [];
yr_k = [];

for in=1:nsim
    in

  % updating arrays
  ur_k = [ur_k, uk_1];
  yr_k = [yr_k, ypk];
  
  % building terms of the OF that depend on current information
  if distModel ~= 0
      el = ysp - Phi*[xmk;de];
  else 
      el = ysp - Phi*xmk;
  end
  
  ct = -el'*Qbar*Theta - uk_1'*Ibar'*Rbar*IM;

  if distModel == 0 || distModel == 1  % solve optimization problem
      % building constraints on the input changes
      Ac_du = [IM;-IM];
      bc_du = [Dumax + Ibar*uk_1;Dumax - Ibar*uk_1];
      
      % solving QP
      [ukk,vkv] = quadprog(H,ct,Ac_du,bc_du,[],[],Umin,Umax);
      V_k = [V_k, vkv];
  else % solve unconstrained problem
        ukk = -H'\ct';
        V_k = [V_k, ukk'*H*ukk + 2*ct*ukk];
  end

  % extracting first input (update uk)
  uk = ukk(1:nu,1);

  % applying input to the plant
  xpk = Ap*xpk + Bp*uk;
  ypk = Cp*xpk;

  % correcting model state with state KF
  % 1. applying input to the model
  xmk = A*xmk + B*uk;
  ymk = C*xmk;

  % 2. compute prediction error
  de = ypk - ymk;

  if distModel == 0
      % 3. correct current model state
      xmk = xmk + Kf*(de);
  end
  
    % loop
    uk_1 = uk;

end













