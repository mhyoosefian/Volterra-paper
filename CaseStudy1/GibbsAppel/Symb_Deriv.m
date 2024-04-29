clc;
clear all; %#ok<*CLALL>
close all;

syms q1 q2 q3 real     % theta1, theta2, x
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms U1 U2 dU1 dU2 real
syms t
syms m1 l m2 g tau P real

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];
U = [U1; U2];
dU = [dU1; dU2];
%% Angular Velocities
omega1 = [0;0;dq1];
omega2 = [0;0;dq2];

%% Transitional Velocities
VA = [dq3; 0; 0];
VG1 = VA + [l/2*dq1*cos(q1); l/2*dq1*sin(q1); 0];
VB = VA + [l*dq1*cos(q1); l*dq1*sin(q1); 0];
VG2 = VB + [l/2*dq2*cos(q2); l/2*dq2*sin(q2); 0];
VC = VB + [l*dq2*cos(q2); l*dq2*sin(q2); 0];

Ptot = m1*VA + m2*VG1 + m2*VG2;
%% Constraints
e = [cos(q2); sin(q2); 0];
f = simplify(dot((VC-VA), e));
a = jacobian(f, dq);
b = 0;

%% Quasi-velocities
u1 = dq1 - dq2;
u2 = dq3;
u = [u1; u2];
Y = jacobian(u, dq);
Z = simplify(u - Y*dq);
W = simplify(inv([Y; a])*([eye(2,2); zeros(1,2)]));
X = -inv([Y; a])*[Z; b];
dq_itu = W*U+X;

%% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VA = simplify(subs(VA, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
VG1 = simplify(subs(VG1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
VG2 = simplify(subs(VG2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
omega1 = simplify(subs(omega1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
omega2 = simplify(subs(omega2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

%% Kinetic Energy
T1 = 0.5*m1*(VA.'*VA);
T2 = 0.5*m2*(VG1.'*VG1) + 0.5*(m2*l*l/12)*(omega1.'*omega1);
T3 = 0.5*m2*(VG2.'*VG2) + 0.5*(m2*l*l/12)*(omega2.'*omega2);
T = simplify(T1 + T2 + T3);

%% Linear accelerations
accelA = zeros(3,1);
for i=1:length(q)
accelA = accelA + diff(VA,q(i))*dq(i);
end
for i=1:length(U)
accelA = accelA + diff(VA,U(i))*dU(i);
end
accelA = (subs(accelA, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

accelG1 = zeros(3,1);
for i=1:length(q)
accelG1 = accelG1 + diff(VG1,q(i))*dq(i);
end
for i=1:length(U)
accelG1 = accelG1 + diff(VG1,U(i))*dU(i);
end
accelG1 = (subs(accelG1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

accelG2 = zeros(3,1);
for i=1:length(q)
accelG2 = accelG2 + diff(VG2,q(i))*dq(i);
end
for i=1:length(U)
accelG2 = accelG2 + diff(VG2,U(i))*dU(i);
end
accelG2 = (subs(accelG2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

%% Angular accelerations
alpha1 = zeros(3,1);
for i=1:length(U)
alpha1 = alpha1 + diff(omega1,U(i))*dU(i);
end
for i=1:length(q)
alpha1 = alpha1 + diff(omega1,q(i))*dq_itu(i);
end
alpha1 = simplify(alpha1);
% alpha1 = (subs(alpha1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

alpha2 = zeros(3,1);
for i=1:length(U)
alpha2 = alpha2 + diff(omega2,U(i))*dU(i);
end
for i=1:length(q)
alpha2 = alpha2 + diff(omega2,q(i))*dq_itu(i);
end
alpha2 = simplify(alpha2);
% alpha2 = (subs(alpha2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

%% Potential Energy
V = -m2*g*l/2*cos(q1) - m2*g*(l*cos(q1) + l/2*cos(q2)) - m1*g*l;

%% Total Mechanical Energy
E = simplify(T + V);

%% Generalized Forces
for i=1:length(U)
Qpr(i,:) = dot(jacobian(omega1, U(i)), [0;0;tau]) + dot(jacobian(omega2, U(i)), [0;0;-tau]) + dot(jacobian(VG1, U(i)), [0;-m2*g;0])...
      + dot(jacobian(VG2, U(i)), [0;-m2*g;0]) + dot(jacobian(VA, U(i)), [0;-m1*g;0]);
end

%% Building Gibss Function
I = m2*l*l/12;
S = 0.5*(m1*(dot(accelA,accelA)) + m2*(dot(accelG1,accelG1) + dot(accelG2,accelG2)))...
    + 0.5*(dot(alpha1,I*alpha1) + dot(alpha2,I*alpha2))...
    + dot(alpha1,cross(omega1, I*omega1)) + dot(alpha2,cross(omega2, I*omega2));

%% Equations of Motion
Sdot1 = diff(S,dU1);
Sdot2 = diff(S,dU2);

eq = [Sdot1;Sdot2] - Qpr;
M_g = simplify(jacobian(eq, dU));
B_g = (eq - M_g*dU);

Ptot_u = m1*VA + m2*VG1 + m2*VG2;
%% Functions
matlabFunction(M_g, 'File', 'M_gibss');
matlabFunction(B_g, 'File', 'B_gibss');
matlabFunction(f, 'File', 'constraints');
matlabFunction(E, 'File', 'Energy_Gibbs');
matlabFunction(Ptot_u, 'File', 'ptot_Gibbs_u');
matlabFunction(Ptot, 'File', 'ptot_Gibbs');
matlabFunction(dq_itu, 'File', 'dq_itu_Gibbs');