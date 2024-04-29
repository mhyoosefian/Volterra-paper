clc;
clear all; %#ok<*CLALL>
close all;

syms q1 q2 q3 real     % theta1, theta2, x
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms t
syms m1 l m2 g tau P real

q = [q1;q2;q3];
dq = [dq1;dq2;dq3];

% Angular Velocities
omega1 = [0;0;dq1];
omega2 = [0;0;dq2];

% Transitional Velocities
VA = [dq3; 0; 0];
VG1 = VA + [l/2*dq1*cos(q1); l/2*dq1*sin(q1); 0];
VB = VA + [l*dq1*cos(q1); l*dq1*sin(q1); 0];
VG2 = VB + [l/2*dq2*cos(q2); l/2*dq2*sin(q2); 0];
VC = VB + [l*dq2*cos(q2); l*dq2*sin(q2); 0];

% Constraint
e = [cos(q2); sin(q2); 0];
f = simplify(dot(VC-VA, e));
a = jacobian(f, dq);
adot = zeros(size(a));
for i =1:length(dq)
   adot = adot + diff(a, q(i))*dq(i); 
end
b = 0;
bdot = 0;

% Kinetic Energy
T1 = 0.5*m1*(VA.'*VA);
T2 = 0.5*m2*(VG1.'*VG1) + 0.5*(m2*l*l/12)*(omega1.'*omega1);
T3 = 0.5*m2*(VG2.'*VG2) + 0.5*(m2*l*l/12)*(omega2.'*omega2);

T = simplify(T1 + T2 + T3);

% Potential Energy
V = -m2*g*l/2*cos(q1) - m2*g*(l*cos(q1) + l/2*cos(q2)) - m1*g*l;
L = T - V;

% Total Mechanical Energy
E = simplify(T + V);

% Generalized Forces
Q = [tau; -tau; 0];

% Lag_IM Equations
pl_pdq = jacobian(L,dq).';
M = jacobian(pl_pdq,dq);
M = simplify(M);

pl_pq = jacobian(L,q).';
B = diff(pl_pdq,t) + jacobian(pl_pdq,q)*dq - pl_pq;
B = simplify(B);
F = Q - B;

A_IM = [eye(3,3) zeros(3,3) zeros(3,1); zeros(3,3) M -a.'; zeros(1,3) -a zeros(1,1)];
B_IM = [dq; F; bdot + adot*dq];
A_IM = simplify(A_IM);
B_IM = simplify(B_IM);

%Linear Momentum
Ptot = m1*VA + m2*(VG1 + VG2);

matlabFunction(A_IM, 'File', 'A');
matlabFunction(B_IM, 'File', 'B');
matlabFunction(f, 'File', 'constraints');
matlabFunction(E, 'File', 'Energy');
matlabFunction(Ptot, 'File', 'Lmom');