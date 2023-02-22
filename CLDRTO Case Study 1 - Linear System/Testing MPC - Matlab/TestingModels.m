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

% plant model
Ap=A;

%changing B matrix: test influence process gain
Bp{1}=B;  
Bp{1}(1)=0.05;  

Bp{2}=B;  
Bp{2}(1)=0.0625; 

Bp{3}=B;  
Bp{3}(1)=0.07;  

Cp=C;

%% Starting simulation
% model and plant state (initial state is known)
xpk = x0;

% model and plant outputs 
ypk = Cp*xpk;

% initial input (for computing input movement du)
uk = 1;

for kk = 1:3
    yr_k{kk} = [];
end

for kk = 1:3
    for in=1:nsim    
        in
      % updating arrays
      yr_k{kk} = [yr_k{kk}, ypk];
      
      % applying input to the plant
      xpk = Ap*xpk + Bp{kk}*uk;
      ypk = Cp*xpk;

    end
end


%% Plotting
figure(1)

    hold on
    plot(yr_k{1},'g:','linewidth',1.5)
    plot(yr_k{2},'b-.','linewidth',1.5)
    plot(yr_k{3},'r','linewidth',1.5)
    
    grid on
    box on
    ylim([0, 1.5])
    xlim([0, nsim])

    ylabel('y')
