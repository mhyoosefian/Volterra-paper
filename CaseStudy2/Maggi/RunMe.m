clc;
clear all %#ok<*CLALL>
simulation_time = [0:0.1:50]; %#ok<*NBRAK>

%% Problem initializing
global m_s m_p a b c Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 tau_p1 tau_p2
m_s = 100; Ixx = 67; Iyy = 67; Izz = 67; Ixy = 5; Ixz = 2; Iyz = 0; a = 2; b = 2; c = 2;
m_p = 10;
tau1 = 0; tau2 = 0; tau3 = 0;  tau_p1 = 0; tau_p2 = 0; 
%% Solving the system using ODE45
options = odeset('maxstep',0.1);
z0 = [pi/10;pi/6;0.03;  0;0;  3;3;9;      0.3;0;-0.3;  0;0.2;  0.1;-0.3;-0.4]; % Main initial condition       %Z=[q;dq]
tic;
[t,z] = ode45(@Dyn_Maggie, simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7); q8 = z(:,8);
dq1 = z(:,9); dq2 = z(:,10); dq3 = z(:,11); dq4 = z(:,12); dq5 = z(:,13); dq6 = z(:,14); dq7 = z(:,15); dq8 = z(:,16);
%% Energy conservation criteria
E = Energy_Maggie(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,b,c,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,m_p,m_s,q1,q2,q3,q4,q5);
E_Err_Mag = (E-E(1))/E(1);
% save('E_Err_Mag');
X = ['Norm of Energy Error = ',num2str(norm((E-E(1))/E(1)))];
disp(X);
%% Linear momentum criteria
LM = Lmom_Maggie(b,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,m_p,m_s,q1,q2,q3,q4,q5);
LMX = LM(1:length(t));
LMY = LM(length(t)+1:2*length(t));
LMZ = LM(2*length(t)+1:3*length(t));
LM_Err_Mag = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Mag');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);