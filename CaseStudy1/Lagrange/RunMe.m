clc; 
close all; 
clear all; %#ok<*CLALL>
global m1 l m2 g tau

m1 = 1; m2 = 0.5;
l = 0.2;
tau = 0;
g = 9.81;
initcond = [pi/2; pi/2; 4; 1; -1; 3];     %[q; dq]  
z0 = [pi/2;pi/2;4;   1;-1;3;  0];     %[q; dq; mu]
simulation_time = 50; %#ok<*NBRAK>
option = odeset('maxstep',0.01);
tic;
[t,z] = ode45(@Dyn_LagIM, [0:0.01:simulation_time],z0,option); 
toc

q1 = z(:,1); q2 = z(:,2); q3 = z(:,3);
dq1 = z(:,4); dq2 = z(:,5); dq3 = z(:,6);



E = Energy(dq1,dq2,dq3,g,l,m1,m2,q1,q2);
E_Err_Lag = (E-E(1))/E(1);
% save('E_Err_Lag ');
X = ['Norm of Energy Error = ',num2str(norm((E-E(1))/E(1)))];
disp(X);

f = constraints(dq1,dq2,l,q1,q2);
f_Err_Lag = f;
% save('f_Err_Lag ');
X = ['Norm of Constraint Error = ',num2str(norm(f))];
disp(X);

% Linear Momentum Conservation
LM = Lmom(dq1,dq2,dq3,l,m1,m2,q1,q2);
LMX = LM(1:length(t),:);
LM_Err_Lag = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Lag ');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);