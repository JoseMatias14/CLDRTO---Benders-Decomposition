% Applying MPC to a linear sytstem with plant-model mismatch
% - MPC without disturbance model
% - MPC with disturbance model (no offset)
% - LQR

% System from paper: Dynamic Real-time Optimization with Closed-Loop
% Prediction (2017) - MZ Jamaludin
clear
close all

%% System configuration 
% Sampling time
T = 1; 
% Simulation time in sampling periods
nsim = 75;
% Number of manipulated inputs
nu = 1;
% Number of controlled outputs
ny = 1;
% Initial value of the inputs
u0= 0.0;
% Initial values of the outputs
y0= 0.0;

%% System model
% Transfer functions
G = tf([1.0],[125.0, 75.0, 15.0, 1.0]);
Gd = c2d(G,T,'zoh');

% state-space model (used in controller)
[A,B,C,D]=ssdata(Gd);
nx=size(A,1);
x0 = zeros(nx,1);

%adding disturbance model
Ad = [A zeros(nx,ny); zeros(ny,nx), eye(ny)];
Bd = [B; zeros(ny,nu)];
Cd = [C, eye(ny)];

% plant model
Ap=A;
Bp=1.1*B;  %multiplying the process matrix by 0.1 to create plant-model mismatch
Cp=C;

%% Controller tuning (see Table 2 - Cases 1a/1b - base case tuning)
% Output prediction horizon
p = 30; 
% Input control horizon 
m = 3;
% Output weights (
q = 1;
% Input weights
r = 1;
% Set-point of the outputs
ys = 1;
% inputs bounds -- no active constraints
umax = 12; %1.2
umin = 0;
dumax = 3; %0.3
%umax = 1.2; %1.2
%umin = 0;
%dumax = 0.3; %0.3

%%%%%%%%%%%%%%%%%%%%
% MPC regular form %
%%%%%%%%%%%%%%%%%%%%
distModel = 0;
[ur1,yr1]=MPC_COMP(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ad,Bd,Cd,Ap,Bp,Cp,umax,umin,dumax,x0,y0,ys,distModel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC with disturbance model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distModel = 1;
[ur2,yr2]=MPC_COMP(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ad,Bd,Cd,Ap,Bp,Cp,umax,umin,dumax,x0,y0,ys,distModel);

%%%%%%%%%%%%%%%%%%%%%
% Unconstrained MPC %
%%%%%%%%%%%%%%%%%%%%%
distModel = 2;
[ur3,yr3]=MPC_COMP(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ad,Bd,Cd,Ap,Bp,Cp,umax,umin,dumax,x0,y0,ys,distModel);

%
%% Plotting
figure(1)

subplot(2,1,1)
    yline(ys,'k--')
    hold on
    plot(yr1(1,1:nsim),'g:','linewidth',1.5)
    plot(yr2(1,1:nsim),'b-.','linewidth',1.5)
    plot(yr3(1,1:nsim),'r','linewidth',1.5)
    
    legend({'SP','MPC','MPC + dist','LQR'},'Location','best')
    grid on
    box on
    ylim([0, 1.5])
    xlim([0, nsim])

    ylabel('y')

subplot(2,1,2)
    hold on
    stairs(ur1(1,1:nsim),'g:','linewidth',1.5)
    stairs(ur2(1,1:nsim),'b-.','linewidth',1.5)
    stairs(ur3(1,1:nsim),'r','linewidth',1.5)
    yline(umax,'k:')
    yline(umin,'k:')
    grid on
    box on
    ylim([-0.1, 2])
    xlim([0, nsim])

    ylabel('u')    
    xlabel('time')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Linear Quadratic Regulator % -----> not used here!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Computing LQR gain 
% R_LQR = 1;
% Q_LQR = ones(4);
% [K_LQR,S_LQR,P_LQR] = lqi(ss(A,B,C,D),Q_LQR,R_LQR);
% 
% % State observer
% Kf = FKalman(ny,A,C,100);
% 
% % starting simulation
% % model and plant state (initial state is known)
% xmk = x0;
% xpk = xmk;
% 
% % computing reference values
% xmr = x0;
% x_temp = x0;
% der = 0;
% umr = 0.9091;
% for ii = 1:100
%     xmr = A*xmr + B*umr; 
%     x_temp = Ap*x_temp + B*umr;
%     der = Cp*x_temp - C*xmr;
% end
% 
% % plant outputs 
% ypk = Cp*xpk;
% ymk = C*xmk;
% de = ypk - ymk;
% 
% % initial input
% uk_1 = u0;
% 
% % for saving values
% ur3 = [];
% yr3 = [];
% 
% for in=1:nsim
%     in
% 
%   % updating arrays
%   ur3 = [ur3, uk_1];
%   yr3 = [yr3, ypk];
%   
%   % using LQR to compute new input
%   uk = -K_LQR*([xmk;de] - [xmr;der]) + umr;
% 
%   % applying input to the plant
%   xpk = Ap*xpk + Bp*uk;
%   ypk = Cp*xpk;
% 
%   % correcting model state with state KF
%   % applying input to the model
%   xmk = A*xmk + B*uk;
%   ymk = C*xmk;
%   
%   % compute prediction error
%   de = ypk - ymk;
% 
%   % loop 
%   uk_1 = uk;
% 
% end
