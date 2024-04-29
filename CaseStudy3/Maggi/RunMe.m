clc;
clear all %#ok<*CLALL>
simulation_time = [0:0.1:50]; %#ok<*NBRAK>

%% Problem initializing
global m_sat l Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 a m
m_sat = 2000; Ixx = 1400; Iyy = 900; Izz = 1100; Ixy = 5; Ixz = 8; Iyz = 3; l = 2.5;
m = 1; a = 0.5;
tau1 = 0; tau2 = 0; tau3 = 0;
pow_0 = Power(0.012,-0.1,0.05,-0.05,0,pi/12,0.08,tau1,tau2,tau3);
z0 = [pi/8;pi/12;0.08; a; 0;0;0;  -0.1;0.05;-0.05; 0; 2;1;0;    pow_0];

%% Solving the system using ODE45
options = odeset('maxstep',0.1);
tic;
[t,z] = ode45(@Dyn_Maggie, simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7);
dq1 = z(:,8); dq2 = z(:,9); dq3 = z(:,10); dq4 = z(:,11); dq5 = z(:,12); dq6 = z(:,13); dq7 = z(:,14);
int_pow = z(:,15);      
%% Energy conservation criteria
F = -0.018*sin(0.089*t) + 0.012*cos(0.0485*t);
E = Energy_Maggie(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,dq1,dq2,dq3,dq4,dq5,dq6,dq7,l,m,m_sat,q1,q2,q3,q4);
Pow_Err_Mag = (E-int_pow);
Pow_Err_Mag = (Pow_Err_Mag - Pow_Err_Mag(1))/Pow_Err_Mag(1);
% save('Pow_Err_Mag');
X = ['Norm of Power Error = ',num2str(norm(Pow_Err_Mag))];
disp(X);

%% Linear momentum criteria
LM = Lmom_Maggie(a,dq1,dq2,dq3,dq4,dq5,dq6,dq7,l,m,m_sat,q1,q2,q3,q4);
LMX = LM(1:length(t));
LMY = LM(length(t)+1:2*length(t));
LMZ = LM(2*length(t)+1:3*length(t));
LM_Err_Mag = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Mag');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);
