% Examples of Generating Collocation Matrix
clear
clc

% Radau Collocation Points 3th degree polynomial
t0=0;
t1=0.15505;
t2=0.64495;
t3 =1.00000;

syms t
l0=(t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3));
l1=(t-t0)*(t-t2)*(t-t3)/((t1-t0)*(t1-t2)*(t1-t3));
l2=(t-t0)*(t-t1)*(t-t3)/((t2-t0)*(t2-t1)*(t2-t3));
l3=(t-t0)*(t-t1)*(t-t2)/((t3-t0)*(t3-t1)*(t3-t2));

% Coefficients for the Lagrange polynomials
L0 = sym2poly(l0);
L1 = sym2poly(l1);
L2 = sym2poly(l2);
L3 = sym2poly(l3);

% Coefficients for the 1st derivatives
dL0 = polyder(L0);
dL1 = polyder(L1);
dL2 = polyder(L2);
dL3 = polyder(L3);

% Collocation matrix: 1st derivatives evaluated at the collocation points
tau = [t0,t1,t2,t3];
adot3 = zeros(4,4);
adot3(1,:) = polyval(dL0,tau);
adot3(2,:) = polyval(dL1,tau);
adot3(3,:) = polyval(dL2,tau);
adot3(4,:) = polyval(dL3,tau);
adot3

% Radau Collocation Points 2nd degree polynomial
t0=0;
t1=0.33333;
t2=1.00000;

syms t
l0=(t-t1)*(t-t2)/((t0-t1)*(t0-t2));
l1=(t-t0)*(t-t2)/((t1-t0)*(t1-t2));
l2=(t-t0)*(t-t1)/((t2-t0)*(t2-t1));

% Coefficients for the Lagrange polynomials
L0 = sym2poly(l0);
L1 = sym2poly(l1);
L2 = sym2poly(l2);

% Coefficients for the 1st derivatives
dL0 = polyder(L0);
dL1 = polyder(L1);
dL2 = polyder(L2);

% Collocation matrix: 1st derivatives evaluated at the collocation points
tau = [t0,t1,t2];
adot2 = zeros(3,3);
adot2(1,:) = polyval(dL0,tau);
adot2(2,:) = polyval(dL1,tau);
adot2(3,:) = polyval(dL2,tau);
adot2

% Radau Collocation Points 1st degree polynomial (backward euler
t0=0;
t1=1.00000;

syms t
l0=(t-t1)/(t0-t1);
l1=(t-t0)/(t1-t0);

% Coefficients for the Lagrange polynomials
L0 = sym2poly(l0);
L1 = sym2poly(l1);

% Coefficients for the 1st derivatives
dL0 = polyder(L0);
dL1 = polyder(L1);

% Collocation matrix: 1st derivatives evaluated at the collocation points
tau = [t0,t1];
adot1 = zeros(2,2);
adot1(1,:) = polyval(dL0,tau);
adot1(2,:) = polyval(dL1,tau);
adot1



