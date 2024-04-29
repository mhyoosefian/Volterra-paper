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
initial_condition = [pi/10;pi/6;0.03;  0;0;  3;3;9;      0.3;0;-0.3;  0;0.2;  0.1;-0.3;-0.4]; % Main initial condition
u = u_Gibbs(0.3,0,-0.3,0,0.2,0.1,-0.3,-0.4,pi/6,0.03);
z0 = [pi/10;pi/6;0.03;  0;0;  3;3;9;     u];
tic;
[t,z] = ode45(@Dyn_Gibbs, simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7); q8 = z(:,8);
U1 = z(:,9); U2 = z(:,10); U3 = z(:,11); U4 = z(:,12); U5 = z(:,13); U6 = z(:,14); U7 = z(:,15); U8 = z(:,16);
%% Energy conservation criteria
E = Energy_Gibbs(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,U5,U6,U7,U8,a,b,c,m_p,m_s,q1,q2,q3,q4,q5);
E_Err_Gibbs = (E-E(1))/E(1);
% save('E_Err_Gibbs');
X = ['Norm of Energy Error = ',num2str(norm((E-E(1))/E(1)))];
disp(X);
%% Linear momentum criteria
LM = Lmom_Gibbs(U1,U2,U3,U4,U5,U6,U7,U8,b,m_p,m_s,q1,q2,q3,q4,q5);
LMX = LM(1:length(t));
LMY = LM(length(t)+1:2*length(t));
LMZ = LM(2*length(t)+1:3*length(t));
LM_Err_Gibbs = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Gibbs');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);
