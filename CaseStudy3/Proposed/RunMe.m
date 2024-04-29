clc;
clear all %#ok<*CLALL>
simulation_time = [0:0.1:50]; %#ok<*NBRAK>

%% Problem initializing
global m_sat l Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 a m px py pz
m_sat = 2000; Ixx = 1400; Iyy = 900; Izz = 1100; Ixy = 5; Ixz = 8; Iyz = 3; l = 2.5;
m = 1; a = 0.5;
tau1 = 0; tau2 = 0; tau3 = 0;

%% Solving the system using ODE45
options = odeset('maxstep',0.1);
initial_condition = [pi/8;pi/12;0.08; a; 0;0;0;  -0.1;0.05;-0.05; 0; 2;1;0]; % Main initial condition
PXYZ = PXYZ_Voltra(a,-0.1,0.05,-0.05,0,2,1,0,l,m,m_sat,pi/8,pi/12,0.08,a);
uNI = U_NI(-0.1,0.05,-0.05,0,pi/12,0.08);
pow_0 = Power(0.012,uNI(1),uNI(2),uNI(3),uNI(4),pi/12,0.08,tau1,tau2,tau3);
px = PXYZ(1); py = PXYZ(2); pz = PXYZ(3);
z0 = [pi/8;pi/12;0.08; a; 0;0;0;     uNI;    pow_0];
tic;
[t,z] = ode45(@Sat_Dynamics_Voltra_IC, simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7);
U1 = z(:,8); U2 = z(:,9); U3 = z(:,10); U4 = z(:,11);
int_pow = z(:,12);    
%% Energy conservation criteria
F = -0.018*sin(0.089*t) + 0.012*cos(0.0485*t);
E = Energy(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,a,l,m,m_sat,px,py,pz,q1,q2,q3,q4);
Pow_Err_Vol = (E-int_pow);
Pow_Err_Vol = (Pow_Err_Vol - Pow_Err_Vol(1))/Pow_Err_Vol(1);
% save('Pow_Err_Vol');
X = ['Norm of Power Error = ',num2str(norm(Pow_Err_Vol))];
disp(X);
%% Linear momentum criteria
LM = Lmoment(U1,U2,U4,a,l,m,m_sat,px,py,pz,q1,q2,q3,q4);
LMX = LM(1)*ones(length(t));
LMY = LM(2)*ones(length(t));
LMZ = LM(3)*ones(length(t));
LM_Err_Vol_IC = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Vol_IC');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);