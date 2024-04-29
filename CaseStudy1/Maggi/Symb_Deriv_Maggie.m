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
ddq = [ddq1;ddq2;ddq3];

% Angular Velocities
omega1 = [0;0;dq1];
omega2 = [0;0;dq2];

% Transitional Velocities
VA = [dq3; 0; 0];
VG1 = VA + [l/2*dq1*cos(q1); l/2*dq1*sin(q1); 0];
VB = VA + [l*dq1*cos(q1); l*dq1*sin(q1); 0];
VG2 = VB + [l/2*dq2*cos(q2); l/2*dq2*sin(q2); 0];
VC = VB + [l*dq2*cos(q2); l*dq2*sin(q2); 0];

% Constraints
e = [cos(q2); sin(q2); 0];
f = simplify(dot((VC-VA), e));
a = jacobian(f, dq);
adot = zeros(size(a));
for i =1:length(dq)
   adot = adot + diff(a, q(i))*dq(i); 
end
b = 0;

% Kinetic Energy
T1 = 0.5*m1*(VA.'*VA);
T2 = 0.5*m2*(VG1.'*VG1) + 0.5*(m2*l*l/12)*(omega1.'*omega1);
T3 = 0.5*m2*(VG2.'*VG2) + 0.5*(m2*l*l/12)*(omega2.'*omega2);

T = T1 + T2 + T3;

% Potential Energy
V = -m2*g*l/2*cos(q1) - m2*g*(l*cos(q1) + l/2*cos(q2)) - m1*g*l;

% Total Mechanical Energy
E = simplify(T + V);

% W Matrix and Quasi-velocities
u1 = dq1 - dq2;
u2 = dq3;
u = [u1; u2];
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = ([Y; a])\([eye(2,2); zeros(1,2)]);
X = -([Y; a])\[Z; b];


% Generalized Forces
for i=1:length(q)
Q(i,:) = dot(jacobian(omega1, dq(i)), [0;0;tau]) + dot(jacobian(omega2, dq(i)), [0;0;-tau]) + dot(jacobian(VG1, dq(i)), [0;-m2*g;0])...
      + dot(jacobian(VG2, dq(i)), [0;-m2*g;0]) + dot(jacobian(VA, dq(i)), [0;-m1*g;0]);
end
Qpr = simplify(W.'*Q);

% Maggies' Equations
pT_pq = jacobian(T, q).';
pT_pdq = jacobian(T, dq).';
d_dt_pT_pdq = zeros(size(dq));
for i=1:length(q)
   d_dt_pT_pdq =  d_dt_pT_pdq + diff(pT_pdq, dq(i))*ddq(i) + diff(pT_pdq, q(i))*dq(i);
end
eq = (W.'*(d_dt_pT_pdq - pT_pq) - Qpr);
M_ = jacobian(eq, ddq);
A_ = [M_; a];
B = simplify(eq - M_*ddq);
B_ = [B; adot*dq];
Ptot = m1*VA + m2*VG1 + m2*VG2;

matlabFunction(A_, 'File', 'A_Maggie');
matlabFunction(B_, 'File', 'B_Maggie');
matlabFunction(f, 'File', 'constraints');
matlabFunction(E, 'File', 'Energy_Maggie');
matlabFunction(Ptot, 'File', 'Lmom_Maggie');
