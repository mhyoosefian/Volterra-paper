clc;
clear all; %#ok<*CLALL>
close all;

syms q1 q2 q3 real     % theta1, theta2, x
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms U dU real
syms t
syms m1 l m2 g tau P real

q = [q1;q2;q3];
q_NI = [q1; q2];
dq = [dq1;dq2;dq3];
dq_NI = [dq1; dq2];

% Angular Velocities
omega1 = [0;0;dq1];
omega2 = [0;0;dq2];

% Transitional Velocities
VA = [dq3; 0; 0];
VG1 = VA + [l/2*dq1*cos(q1); l/2*dq1*sin(q1); 0];
VB = VA + [l*dq1*cos(q1); l*dq1*sin(q1); 0];
VG2 = VB + [l/2*dq2*cos(q2); l/2*dq2*sin(q2); 0];
VC = VB + [l*dq2*cos(q2); l*dq2*sin(q2); 0];
PXYZ = simplify(m1*VA + m2*(VG1 + VG2));

% Kinetic Energy
T1 = 0.5*m1*(VA.'*VA);
T2 = 0.5*m2*(VG1.'*VG1) + 0.5*(m2*l*l/12)*(omega1.'*omega1);
T3 = 0.5*m2*(VG2.'*VG2) + 0.5*(m2*l*l/12)*(omega2.'*omega2);
T = simplify(T1 + T2 + T3);

% Potential Energy
V = -m2*g*l/2*cos(q1) - m2*g*(l*cos(q1) + l/2*cos(q2)) - m1*g*l;
% E = simplify(T+V);
L = T - V;

% Linear Momentum
% Ptot = m1*VA + m2*(VG1 + VG2);

% Constraints
e = [cos(q2); sin(q2); 0];
f = simplify(dot((VC-VA), e));
a = jacobian(f, dq);
b = 0;

% W Matrix and Quasi-velocities
pl_pdq = jacobian(L,dq).';
M = jacobian(pl_pdq,dq);

M11 = M(1:2, 1:2);
M12 = M(1:2,3);
M21 = M12.';
M22 = M(3,3);
M2 = simplify([M21 M22]);
N2 = simplify(diff(L, dq3) - M21*dq_NI - M22*dq3);

u = dq2 - dq1;
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = simplify(([Y; M2; a])\([eye(1,1); zeros(1,1); zeros(1,1)]));
X = simplify(-([Y; M2; a])\[Z;N2-P ;b]);
dq_itu = simplify(W*U+X);

% for i=1:length(q)
% Q(i,:) = dot(jacobian(omega1, dq(i)), [0;0;tau]) + dot(jacobian(omega2, dq(i)), [0;0;-tau]) + dot(jacobian(VG1, dq(i)), [0;-m2*g;0])...
%       + dot(jacobian(VG2, dq(i)), [0;-m2*g;0]) + dot(jacobian(VA, dq(i)), [0;-m1*g;0]);
% end
% Qpr2 = simplify(W.'*Q);
% Qpr2 = simplify(subs(Qpr2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VA = (subs(VA, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
VG1 = (subs(VG1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
VB = (subs(VB, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
VG2 = (subs(VG2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
omega1 = (subs(omega1, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
omega2 = (subs(omega2, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

% Kinetic Energy in terms of U
T1 = 0.5*m1*(VA.'*VA);
T2 = 0.5*m2*(VG1.'*VG1) + 0.5*(m2*l*l/12)*(omega1.'*omega1);
T3 = 0.5*m2*(VG2.'*VG2) + 0.5*(m2*l*l/12)*(omega2.'*omega2);
T = simplify(T1 + T2 + T3);
E = simplify(T + V);
% Generalized Forces
% Q = [-tau - 3*m2*g*l*sin(q1)/2; tau - m2*g*l*sin(q2)/2; 0];
% Qpr = (W.'*Q);
Qpr = dot(jacobian(omega1, U), [0;0;tau]) + dot(jacobian(omega2, U), [0;0;-tau]) + dot(jacobian(VG1, U), [0;-m2*g;0])...
      + dot(jacobian(VG2, U), [0;-m2*g;0]) + dot(jacobian(VA, U), [0;-m1*g;0]);

% Reduced Equations of Voltra
pT_pU = simplify(jacobian(T, U).');
d_dt_pT_pU = zeros(size(pT_pU));
for i =1:length(q)
    d_dt_pT_pU = d_dt_pT_pU + diff(pT_pU, q(i))*dq_itu(i);
end
d_dt_pT_pU = simplify(d_dt_pT_pU  + diff(pT_pU, U)*dU);
d_dt_pT_pU = simplify(subs(d_dt_pT_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

pVA_pU = jacobian(VA, U);
pVG1_pU = jacobian(VG1, U);
pVG2_pU = jacobian(VG2, U);
d_dt_pVA_pU = zeros(size(pVA_pU));
d_dt_pVG1_pU = zeros(size(pVG1_pU));
d_dt_pVG2_pU = zeros(size(pVG2_pU));
for i =1:length(q)
    d_dt_pVA_pU = d_dt_pVA_pU + diff(pVA_pU, q(i))*dq_itu(i);
    d_dt_pVG1_pU = d_dt_pVG1_pU + diff(pVG1_pU, q(i))*dq_itu(i);
    d_dt_pVG2_pU = d_dt_pVG2_pU + diff(pVG2_pU, q(i))*dq_itu(i);
end
d_dt_pVA_pU = simplify(d_dt_pVA_pU  + diff(pVA_pU, U)*dU);
d_dt_pVG1_pU = simplify(d_dt_pVG1_pU  + diff(pVG1_pU, U)*dU);
d_dt_pVG2_pU = simplify(d_dt_pVG2_pU  + diff(pVG2_pU, U)*dU);
d_dt_pVA_pU = simplify(subs(d_dt_pVA_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
d_dt_pVG1_pU = simplify(subs(d_dt_pVG1_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
d_dt_pVG2_pU = simplify(subs(d_dt_pVG2_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

P1 = m1*VA;
P2 = m2*VG1;
P3 = m2*VG2;


pomega1_pU = jacobian(omega1, U);
pomega2_pU = jacobian(omega2, U);
d_dt_pomega1_pU = zeros(size(pomega1_pU ));
d_dt_pomega2_pU = zeros(size(pomega2_pU ));
for i =1:length(q)
    d_dt_pomega1_pU = d_dt_pomega1_pU + diff(pomega1_pU, q(i))*dq_itu(i);
    d_dt_pomega2_pU = d_dt_pomega2_pU + diff(pomega2_pU, q(i))*dq_itu(i);
end
d_dt_pomega1_pU = simplify(d_dt_pomega1_pU  + diff(pomega1_pU, U)*dU);
d_dt_pomega2_pU = simplify(d_dt_pomega2_pU  + diff(pomega2_pU, U)*dU);
d_dt_pomega1_pU = simplify(subs(d_dt_pomega1_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));
d_dt_pomega2_pU = simplify(subs(d_dt_pomega2_pU, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

H2 = 1/12*m2*l*l*omega1;
H3 = 1/12*m2*l*l*omega2;


% eq = (d_dt_pT_pU - simplify(dot(P1, d_dt_pVA_pU)) - simplify(dot(P2, d_dt_pVG1_pU)) - simplify(dot(P3, d_dt_pVG2_pU))...
% - simplify(dot(H2, d_dt_pomega1_pU)) - simplify(dot(H3, d_dt_pomega2_pU)) - Qpr);
eq = (d_dt_pT_pU - simplify(d_dt_pVA_pU.'*P1) - simplify(d_dt_pVG1_pU.'*P2) - simplify(d_dt_pVG2_pU.'*P3)...
- simplify(d_dt_pomega1_pU.'*H2) - simplify(d_dt_pomega2_pU.'*H3) - Qpr);
A_ = simplify(jacobian(eq, dU));
B_ = simplify((eq - A_*dU));


% Linear Momentum
Ptot = m1*VA + m2*(VG1 + VG2);

% Potential Energy
f = (subs(f, {dq1 dq2 dq3}, {dq_itu(1) dq_itu(2) dq_itu(3)}));

% Save
matlabFunction(A_, 'File', 'A_Voltra');
matlabFunction(B_, 'File', 'B_Voltra');
% matlabFunction(AI, 'File', 'AI_Voltra');
% matlabFunction(bI, 'File', 'BI_Voltra');
matlabFunction(f, 'File', 'constraints');
matlabFunction(E, 'File', 'Energy_Voltra');
matlabFunction(Ptot, 'File', 'Lmom_Voltra');
% matlabFunction(W, 'File', 'W_Voltra');
matlabFunction(dq_itu, 'File', 'dq_itu_Voltra');
matlabFunction(M21, 'File', 'M21_Voltra');
matlabFunction(M22, 'File', 'M22_Voltra');
matlabFunction(PXYZ,'File','PXYZ_Voltra');