clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi  gamma1   gamma2    X     Y     Z  
syms q1     q2      q3    q4        q5      q6    q7    q8    real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 real
syms m_s m_p Ixx Iyy Izz Ixy Ixz Iyz a b c real
syms tau1 tau2 tau3 tau_p1 tau_p2 real
syms U1 U2 U3 U4 U5 real
syms dU1 dU2 dU3 dU4 dU5 real
syms px py pz real

P = [px; py; pz];

U = [U1; U2; U3; U4; U5];
dU = [dU1; dU2; dU3; dU4; dU5];
q = [q1;q2;q3;q4;q5;q6;q7;q8];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8]; 
dq_I = [dq6;dq7;dq8];
dq_NI = [dq1;dq2;dq3;dq4;dq5];

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
VG1 = simplify(r1.'*VG1);
VG1_XYZ = simplify(R.'*(VG + cross(omega_s, r_AG) + r1.'*cross(omega_p1,r_G1A)));

r_BG = [b/2; a/2; 0];
VB = simplify(r2*(VG + cross(omega_s, r_BG))); %expressed in x2y2z2 (panel2 coordinate system)
r_G2B = [-b/2; 0; 0];
VG2 = simplify(VB + cross(omega_p2,r_G2B)); %expressed in x2y2z2 (panel2 coordinate system)
VG2 = simplify(r2.'*VG2);
VG2_XYZ = simplify(R.'*(VG + cross(omega_s, r_BG) + r2.'*cross(omega_p2,r_G2B)));

PXYZ = simplify(m_s*VG_XYZ + m_p*(VG1_XYZ + VG2_XYZ));
%% Moments of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)
I_p1 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x1y1z1 (panel1 coordinate system)
I_p2 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x2y2z2 (panel2 coordinate system)

%% Kinetic Energy
T_s = simplify(1/2*m_s*VG.'*VG + 1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_p1 = simplify(1/2*m_p*VG1.'*VG1 + 1/2*omega_p1.'*(I_p1*omega_p1));
T_p2 = simplify(1/2*m_p*VG2.'*VG2 + 1/2*omega_p2.'*(I_p2*omega_p2));
T = simplify(T_s + T_p1 + T_p2);

%% Quasi-velocities
pT_pdq = jacobian(T, dq).';
M = (jacobian(pT_pdq, dq));
M11 = M(1:5,1:5);
M12 = M(1:5,6:8);
M21 = M12.';
M22 = M(6:8,6:8);
M2 = simplify([M21 M22]);
N2 = simplify(jacobian(T, dq_I).' - M21*dq_NI - M22*dq_I);

u1 = omega_s(1);
u2 = omega_s(2);
u3 = omega_s(3);
u4 = dq4;
u5 = dq5;
u = [u1; u2; u3; u4; u5];
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = simplify(inv([Y; M2])*[eye(5,5); zeros(3,5)]);
X = simplify(-inv([Y; M2])*[Z; N2-P]);
dq_itu = simplify(W*U+X);

%% Generalized Forces
Q = jacobian(omega_s, dq).'*[tau1; tau2; tau3] + jacobian(omega_p1, dq).'*[0;0;tau_p1] + jacobian(omega_p2, dq).'*[0;0;tau_p2];
Q = simplify(Q);
Qpr = W.'*Q;

%% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VG = simplify(subs(VG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG1 = simplify(subs(VG1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG2 = simplify(subs(VG2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG_XYZ = simplify(subs(VG_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG1_XYZ = simplify(subs(VG1_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
VG2_XYZ = simplify(subs(VG2_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_s = simplify(subs(omega_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_p1 = simplify(subs(omega_p1, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
omega_p2 = simplify(subs(omega_p2, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

%% Kinetic Energy 
T = simplify(subs(T, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

%% Reduced Equations of Voltra
pT_pU = simplify(jacobian(T, U).');
d_dt_pT_puNI = zeros(size(pT_pU));
for i =1:length(q)
    d_dt_pT_puNI = d_dt_pT_puNI + diff(pT_pU, q(i))*dq_itu(i);
end
d_dt_pT_puNI = simplify(d_dt_pT_puNI);
for i =1:length(U)
    d_dt_pT_puNI = d_dt_pT_puNI + diff(pT_pU, U(i))*dU(i);
end
d_dt_pT_puNI = simplify(subs(d_dt_pT_puNI, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

pVG_pU = jacobian(VG, U);
pVG1_pU = jacobian(VG1, U);
pVG2_pU = jacobian(VG2, U);
d_dt_pVG_pU = zeros(size(pVG_pU));
d_dt_pVG1_pU = zeros(size(pVG1_pU));
d_dt_pVG2_pU = zeros(size(pVG2_pU));
for i =1:length(q)
    d_dt_pVG_pU = d_dt_pVG_pU + diff(pVG_pU, q(i))*dq_itu(i);
    d_dt_pVG1_pU = d_dt_pVG1_pU + diff(pVG1_pU, q(i))*dq_itu(i);
    d_dt_pVG2_pU = d_dt_pVG2_pU + diff(pVG2_pU, q(i))*dq_itu(i);
end
d_dt_pVG_pU = simplify(d_dt_pVG_pU);
d_dt_pVG1_pU = simplify(d_dt_pVG1_pU);
d_dt_pVG2_pU = simplify(d_dt_pVG2_pU);
for i =1:length(U)
    d_dt_pVG_pU = d_dt_pVG_pU + diff(pVG_pU, U(i))*dU(i);
    d_dt_pVG1_pU = d_dt_pVG1_pU + diff(pVG1_pU, U(i))*dU(i);
    d_dt_pVG2_pU = d_dt_pVG2_pU + diff(pVG2_pU, U(i))*dU(i);
end
d_dt_pVG_pU = simplify(subs(d_dt_pVG_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
d_dt_pVG1_pU = simplify(subs(d_dt_pVG1_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
d_dt_pVG2_pU = simplify(subs(d_dt_pVG2_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
for i=1:length(U)
    d_dt_pVG_pU(:,i) = simplify(d_dt_pVG_pU(:,i) + cross(omega_s, pVG_pU(:,i)));
    d_dt_pVG1_pU(:,i) = simplify(d_dt_pVG1_pU(:,i) + cross(omega_s, pVG1_pU(:,i)));
    d_dt_pVG2_pU(:,i) = simplify(d_dt_pVG2_pU(:,i) + cross(omega_s, pVG2_pU(:,i)));
end
P1 = m_s*VG;
P2 = m_p*VG1;
P3 = m_p*VG2;


pomegas_pU = jacobian(omega_s, U);
pomegap1_pU = jacobian(omega_p1, U);
pomegap2_pU = jacobian(omega_p2, U);
d_dt_pomegas_pU = zeros(size(pomegas_pU));
d_dt_pomegap1_pU = zeros(size(pomegap1_pU));
d_dt_pomegap2_pU = zeros(size(pomegap2_pU));
for i =1:length(q)
    d_dt_pomegas_pU = d_dt_pomegas_pU + diff(pomegas_pU, q(i))*dq_itu(i);
    d_dt_pomegap1_pU = d_dt_pomegap1_pU + diff(pomegap1_pU, q(i))*dq_itu(i);
    d_dt_pomegap2_pU = d_dt_pomegap2_pU + diff(pomegap2_pU, q(i))*dq_itu(i);
end
d_dt_pomegas_pU = simplify(d_dt_pomegas_pU);
d_dt_pomegap1_pU = simplify(d_dt_pomegap1_pU);
d_dt_pomegap2_pU = simplify(d_dt_pomegap2_pU);
for i =1:length(U)
    d_dt_pomegas_pU = d_dt_pomegas_pU + diff(pomegas_pU, U(i))*dU(i);
    d_dt_pomegap1_pU = d_dt_pomegap1_pU + diff(pomegap1_pU, U(i))*dU(i);
    d_dt_pomegap2_pU = d_dt_pomegap2_pU + diff(pomegap2_pU, U(i))*dU(i);
end
d_dt_pomegas_pU = simplify(subs(d_dt_pomegas_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
d_dt_pomegap1_pU = simplify(subs(d_dt_pomegap1_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));
d_dt_pomegap2_pU = simplify(subs(d_dt_pomegap2_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7) dq_itu(8)}));

for i=1:length(U)
    d_dt_pomegas_pU(:,i) = simplify(d_dt_pomegas_pU(:,i) + cross(omega_s, pomegas_pU(:,i)));
    d_dt_pomegap1_pU(:,i) = simplify(d_dt_pomegap1_pU(:,i) + cross(omega_p1, pomegap1_pU(:,i)));
    d_dt_pomegap2_pU(:,i) = simplify(d_dt_pomegap2_pU(:,i) + cross(omega_p2, pomegap2_pU(:,i)));
end

H1 = I_s*omega_s; %expressed in xyz(satelite coordinate system)
H2 = I_p1*omega_p1; %expressed in x1y1z1 (panel1 coordinate system)
H3 = I_p2*omega_p2; %expressed in x2y2z2 (panel2 coordinate system)

eq = (d_dt_pT_puNI - simplify(d_dt_pVG_pU.'*P1) - simplify(d_dt_pVG1_pU.'*P2) - simplify(d_dt_pVG2_pU.'*P3)...
- simplify(d_dt_pomegas_pU.'*H1) - simplify(d_dt_pomegap1_pU.'*H2) - simplify(d_dt_pomegap2_pU.'*H3) - Qpr);
eq = simplify(eq);
A_ = simplify(jacobian(eq, dU));
B_ = simplify((eq - A_*dU));


%% Linear momentum
Ptot = simplify(m_s*VG_XYZ + m_p*(VG1_XYZ + VG2_XYZ));

%% Functions production
matlabFunction(T,'File','Energy');
matlabFunction(Ptot,'File','Lmoment');
matlabFunction(A_,'File','MassMat');
matlabFunction(B_,'File','Bias');
matlabFunction(dq_itu,'File','dq_itu_Voltra');
matlabFunction([u1; u2; u3; u4; u5],'File','U_NI');     %Is used for initial state vector
matlabFunction(PXYZ,'File','PXYZ_Voltra');
matlabFunction(M21,'File','M21_Voltra');
matlabFunction(M22,'File','M22_Voltra');