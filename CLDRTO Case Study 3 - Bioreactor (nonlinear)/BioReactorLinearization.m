close all
clear 
clc

% for example:
import casadi.*

% parameter values
Vm_val = 0.75; %[1/h]
Ks = 0.5; %[g/L]
D = 0.5; %[1/h] --> D = F/V, Fbar = 2 [L/h], V = 4 [L] 
% F is fixed in the example
h = 1.0; %[h] integration step
time = 0:h:(20*h);

% steady state values
% inputs 
CinBar = 1; % [g/L]

% states
% obtained analytically
Cbar = (CinBar - Ks - Vm_val/D)/2 + sqrt((CinBar - Ks - Vm_val/D)^2 +4*CinBar*Ks)/2;
%Cbar = 0.358; %[g/L] for nominal
Pbar = (1/D)*(Vm_val*Cbar/(Ks + Cbar));
%Pbar = 0.642; %[g/L] for nominal

% inputs for simulation
uSim = 1.5; % 1.1 | 1.5

%% Nonlinear model
% Declare model variables
C = SX.sym('C'); % substrate (reactant)
P = SX.sym('P'); % biomass (product)
states = [C;P];

% Declare system inputs
Cin = SX.sym('Cin'); % inlet concentration

% Declare parameters
Vm_par = SX.sym('Vm_par');
p = [Cin;Vm_par];

%  Monod
    mu = Vm_par*C/(Ks + C);
  
% dynamic equations
    rhs = [D*(Cin - C) - mu;
            mu - D*P];
        
% CVODES
systemStruct = struct('x',states,'p',p,'ode',rhs);
% choosing step
opts = struct('tf',h);
F = integrator('F','cvodes',systemStruct,opts);

%initial condition
xk = [Cbar;Pbar];

%for ploting
xCvodes = xk;

%loop in time
for k = 1:(length(time) - 1)
   %Integrating using cvodes
   I_out = F('x0',xk,'p',[uSim;Vm_val]);
   xk = full(I_out.xf);
   
   xCvodes = [xCvodes, xk];

end

%ploting to compare the results Euler vs Cvodes
figure(1)
plot(time,xCvodes(1,:),'rx',time,xCvodes(2,:),'bx')

%% Obtaining the affine version of the model 
%%%%%%%%%%%%%%%%%%%%
% EVALUATING MODEL %
%%%%%%%%%%%%%%%%%%%%
f_x = Function('f_x',{states,p},{rhs});

% calculating them analytically 
Jx = [-D - Vm_val*Ks/(Ks + Cbar)^2, 0;
    Vm_val*Ks/(Ks + Cbar)^2, -D];

Ju = [D;0];

%%%%%%%%%%%%%%%%%
% SENSITIVITIES %
%%%%%%%%%%%%%%%%%
Sxx = F.factory('sensParStates',{'x0','p'},{'jac:xf:x0'});
Sxu = F.factory('sensParStates',{'x0','p'},{'jac:xf:p'});

Jx_sens = full(Sxx([Cbar;Pbar],[CinBar;Vm_val]));
Ju_sens = full(Sxu([Cbar;Pbar],[CinBar;Vm_val]));

% creating state space model 
sys = ss(Jx,Ju,[0,1],0);

% converting to discrete
d_sys = c2d(sys,h);

%initial condition (deviation from bar)
xkL = [0;0];
xkE = [Cbar;Pbar];
xkS = [Cbar;Pbar];
I_out = F('x0',xkS,'p',[1.0;Vm_val]);
xfbar = full(I_out.xf);

%for ploting
xLin = xkL + [Cbar;Pbar];
xEuler = xkE;
xSens = xkS;

%loop in time
for k = 1:(length(time) - 1)
   %"Integrating" using linearized model
   xkL = d_sys.A*xkL + d_sys.B*(uSim - CinBar);
   xLin = [xLin, xkL + [Cbar;Pbar]]; %  recovering original state

   %"Integrating" using EULER
   xkE = xkE + full(f_x(xkE,[uSim;Vm_val]))*h;
   xEuler = [xEuler, xkE]; 

   % using sensitivities
   xkS = xfbar + Jx_sens*(xkS - [Cbar;Pbar]) + Ju_sens(:,1)*(uSim - CinBar);
   xSens = [xSens, xkS]; 
end

figure(1)
hold on 
plot(time,xLin(1,:),'r-',time,xLin(2,:),'b-')
plot(time,xSens(1,:),'r:',time,xSens(2,:),'b:')
plot(time,xEuler(1,:),'r-.',time,xEuler(2,:),'b-.')
grid on
legend('S_{int}','B_{int}','S_{lin}','B_{lin}','S_{sen}','B_{sen}','S_{eul}','B_{eul}','Location','southeast')
xlabel('t [h]')
ylabel('C,P [g/L]')


%% Building Model Uncertainty
Vm_nom = Vm_val; %[1/h]
Vm_min = Vm_val*0.8;
Vm_max = Vm_val*1.2;

% number of scenarios
nScen = 100;
        
Vm_array = linspace(Vm_min, Vm_max, nScen);

% array for saving models
modelArray = [];

% Jacobian function
% calculating the jacobian w.r.t. states (i.e. derivative of f w.r.t. x)
Jx_AD = jacobian(rhs,states); 
Ju_AD = jacobian(rhs,Cin); 

f_x = Function('f_x',{states,p},{Jx_AD}); 
f_u = Function('f_u',{states,p},{Ju_AD}); 

for jj = 1:nScen

    % computing matrices
    Jx_AD_num = full(f_x([Cbar;Pbar],[CinBar;Vm_array(jj)]));
    Ju_AD_num = full(f_u([Cbar;Pbar],[CinBar;Vm_array(jj)]));

    % creating state space model 
    sys = ss(Jx_AD_num,Ju_AD_num,[0,1],0);
    
    % converting to discrete
    d_sys = c2d(sys,h);

    modelArray = [modelArray; [d_sys.A, d_sys.B]];

end

csvwrite('BioreactorModel',modelArray)
% save('BioreactorModel',"modelArray")

%% Building CSTR in series model 
% testing number of reactors
nReactor = 5;

% array for saving the values
for nn = 1:nReactor
    yLinMult{nn} = [];
end

% prepare plot
figure(2)
plotColor = colormap(parula(nReactor));
% label for plotting
labelMult = {};
for nn = 1:nReactor
    labelMult = {labelMult{:}, ['nr = ', num2str(nn)]};
end

for nn = 1:5 %nReactor
   
    % building multiple CSTR model
    Amult = kron(eye(nn),d_sys.A);
    Bmult = repmat(d_sys.B,[nn,1]);
    Cmult = 1/nn.*repmat([0,1],[1,nn]);
    
    % preparing loop in time
    xkL = repmat([0;0],[nn,1]);

    % saving msasurement for ploting
    yLinMult{nn} = Cmult*xkL + Pbar;

    for k = 1:(length(time) - 1)
       %"Integrating" using linearized model
       xkL = Amult*xkL + Bmult*(uSim - CinBar);
       yLinMult{nn} = [yLinMult{nn}, Cmult*xkL + Pbar]; %  recovering nondeviation measurement

    end
   
    hold on 
    plot(time,yLinMult{nn},'color',plotColor(nn,:),'LineWidth',1.5)
    
end
    
grid on
legend(labelMult,'Location','best')
xlabel('t [h]')
ylabel('P [g/L]')


