clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi  gamma1   gamma2    X     Y     Z  
syms q1     q2      q3    q4        q5      q6    q7    q8    real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 real
syms m_s m_p Ixx Iyy Izz Ixy Ixz Iyz a b c real
syms tau1 tau2 tau3 tau_p1 tau_p2 real
syms U1 U2 U3 U4 U5 U6 U7 U8 real
syms dU1 dU2 dU3 dU4 dU5 dU6 dU7 dU8 real

U = [U1; U2; U3; U4; U5; U6; U7; U8];
dU = [dU1; dU2; dU3; dU4; dU5; dU6; dU7; dU8];
q = [q1;q2;q3;q4;q5;q6;q7;q8];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8]; 

Rphi = [1 0 0;0 cos(q3) sin(q3);0 -sin(q3) cos(q3)];
Rtheta = [cos(q2) 0 -sin(q2);0 1 0;sin(q2) 0 cos(q2)];
Rpsi = [cos(q1) sin(q1) 0;-sin(q1) cos(q1) 0;0 0 1];
R = simplify(Rphi*Rtheta*Rpsi);

r1 = [cos(q4) sin(q4) 0; -sin(q4) cos(q4) 0;0 0 1];
r2 = [cos(q5) -sin(q5) 0; sin(q5) cos(q5) 0;0 0 1];

%% Angular velocities
omega_s = simplify(Rphi*Rtheta*[0;0;dq1] + Rphi*[0;dq2;0] + [dq3;0;0]); %expressed in xyz (satelite coordinate system)
omega_p1 = simplify(r1*(omega_s + [0;0;dq4])); %expressed in x1y1z1 (panel1 coordinate system)
omega_p2 = simplify(r2*(omega_s + [0;0;-dq5])); %expressed in x2y2z2 (panel2 coordinate system)

%% Linear velocities
VG_XYZ = [dq6; dq7; dq8];
VG = simplify(R*VG_XYZ); %expressed in xyz(satelite coordinate system)
r_AG = [b/2; -a/2; 0];
VA = simplify(r1*(VG + cross(omega_s, r_AG))); %expressed in x1y1z1 (panel1 coordinate system)
r_G1A = [-b/2; 0; 0];
VG1 = simplify(VA + cross(omega_p1,r_G1A)); %expressed in x1y1z1 (panel1 coordinate system)
% VG1 = simplify(r1.'*VG1);
VG1_XYZ = simplify(R.'*(VG + cross(omega_s, r_AG) + r1.'*cross(omega_p1,r_G1A)));

r_BG = [b/2; a/2; 0];
VB = simplify(r2*(VG + cross(omega_s, r_BG))); %expressed in x2y2z2 (panel2 coordinate system)
r_G2B = [-b/2; 0; 0];
VG2 = simplify(VB + cross(omega_p2,r_G2B)); %expressed in x2y2z2 (panel2 coordinate system)
% VG2 = simplify(r2.'*VG2);
VG2_XYZ = simplify(R.'*(VG + cross(omega_s, r_BG) + r2.'*cross(omega_p2,r_G2B)));

Ptot= simplify(m_s*VG_XYZ + m_p*(VG1_XYZ + VG2_XYZ));

%% Moments of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)
I_p1 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x1y1z1 (panel1 coordinate system)
I_p2 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x2y2z2 (panel2 coordinate system)

%% Quasi-velocities
u1 = omega_s(1);
u2 = omega_s(2);
u3 = omega_s(3);
u4 = dq4;
u5 = dq5;
u6 = dq6;
u7 = dq7;
u8 = dq8;
u = [u1; u2; u3; u4; u5; u6; u7; u8];
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = simplify(inv(Y));
X = simplify(-Y\Z);
dq_itu = W*U+X;

%% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VG = (subs(VG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG1 = (subs(VG1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG2 = (subs(VG2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG_XYZ = (subs(VG_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG1_XYZ = (subs(VG1_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG2_XYZ = (subs(VG2_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_s = (subs(omega_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_p1 = (subs(omega_p1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_p2 = (subs(omega_p2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

%% Kinetic energy
T_s = simplify(1/2*m_s*VG.'*VG + 1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_p1 = simplify(1/2*m_p*VG1.'*VG1 + 1/2*omega_p1.'*(I_p1*omega_p1));
T_p2 = simplify(1/2*m_p*VG2.'*VG2 + 1/2*omega_p2.'*(I_p2*omega_p2));
T = simplify(T_s + T_p1 + T_p2);

%% Linear accelerations
accelG = zeros(3,1);
for i=1:length(q)
accelG = accelG + diff(VG,q(i))*dq(i);
end
for i=1:length(U)
accelG = accelG + diff(VG,U(i))*dU(i);
end
accelG = simplify(accelG + cross(omega_s, VG));
accelG = (subs(accelG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

accelG1 = zeros(3,1);
for i=1:length(q)
accelG1 = accelG1 + diff(VG1,q(i))*dq(i);
end
for i=1:length(U)
accelG1 = accelG1 + diff(VG1,U(i))*dU(i);
end
accelG1 = simplify(accelG1 + cross(omega_p1, VG1));
accelG1 = (subs(accelG1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

accelG2 = zeros(3,1);
for i=1:length(q)
accelG2 = accelG2 + diff(VG2,q(i))*dq(i);
end
for i=1:length(U)
accelG2 = accelG2 + diff(VG2,U(i))*dU(i);
end
accelG2 = simplify(accelG2 + cross(omega_p2, VG2));
accelG2 = (subs(accelG2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

%% Angular accelerations
alpha_s = zeros(3,1);
for i=1:length(q)
alpha_s = alpha_s + diff(omega_s,q(i))*dq(i);
end
for i=1:length(U)
alpha_s = alpha_s + diff(omega_s,U(i))*dU(i);
end
alpha_s = simplify(alpha_s);
alpha_s = (subs(alpha_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

alpha_p1 = zeros(3,1);
for i=1:length(q)
alpha_p1 = alpha_p1 + diff(omega_p1,q(i))*dq(i);
end
for i=1:length(U)
alpha_p1 = alpha_p1 + diff(omega_p1,U(i))*dU(i);
end
alpha_p1 = simplify(alpha_p1);
alpha_p1 = (subs(alpha_p1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

alpha_p2 = zeros(3,1);
for i=1:length(q)
alpha_p2 = alpha_p2 + diff(omega_p2,q(i))*dq(i);
end
for i=1:length(U)
alpha_p2 = alpha_p2 + diff(omega_p2,U(i))*dU(i);
end
alpha_p2 = simplify(alpha_p2);
alpha_p2 = (subs(alpha_p2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

%% Generalized Forces
Q = jacobian(omega_s, U).'*[tau1; tau2; tau3] + jacobian(omega_p1, U).'*[0;0;tau_p1] + jacobian(omega_p2, U).'*[0;0;tau_p2];
Q = simplify(Q);
Qpr = W.'*Q;

%% Building Gibss Function
S = 0.5*(m_s*(dot(accelG,accelG)) + m_p*(dot(accelG1,accelG1) + dot(accelG2,accelG2)))...
    + 0.5*(dot(alpha_s,I_s*alpha_s) + dot(alpha_p1,I_p1*alpha_p1) + dot(alpha_p2,I_p2*alpha_p2))...
    + dot(alpha_s,cross(omega_s, I_s*omega_s)) + dot(alpha_p1,cross(omega_p1, I_p1*omega_p1))...
    + dot(alpha_p2,cross(omega_p2, I_p2*omega_p2));

%% Equations of Motion
Sdot1 = diff(S,dU1);
Sdot2 = diff(S,dU2);
Sdot3 = diff(S,dU3);
Sdot4 = diff(S,dU4);
Sdot5 = diff(S,dU5);
Sdot6 = diff(S,dU6);
Sdot7 = diff(S,dU7);
Sdot8 = diff(S,dU8);

eq = [Sdot1;Sdot2;Sdot3;Sdot4;Sdot5;Sdot6;Sdot7;Sdot8] - Qpr;
M_g = simplify(jacobian(eq, dU));
B_g = simplify(eq - M_g*dU);

%% Functions
matlabFunction(M_g, 'File', 'M_gibss');
matlabFunction(B_g, 'File', 'B_gibss');
matlabFunction(T, 'File', 'Energy_Gibbs');
matlabFunction(Ptot, 'File', 'ptot_Gibbs');
matlabFunction(dq_itu, 'File', 'dq_itu_Gibbs');
matlabFunction(u, 'File', 'u_Gibbs');     % Wii be used in initial condition