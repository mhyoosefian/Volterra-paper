clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi  gamma1   gamma2    X     Y     Z  
syms q1     q2      q3    q4        q5      q6    q7    q8    real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 real
syms m_s m_p Ixx Iyy Izz Ixy Ixz Iyz a b c real
syms tau1 tau2 tau3 tau_p1 tau_p2 real

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

%% Moments of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)
I_p1 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x1y1z1 (panel1 coordinate system)
I_p2 = (1/12*m_p)*[c^2 0 0;0 b^2+c^2 0;0 0 b^2]; %expressed in x2y2z2 (panel2 coordinate system)

%% Lagrangian function
T_s = simplify(1/2*m_s*VG.'*VG + 1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_p1 = simplify(1/2*m_p*VG1.'*VG1 + 1/2*omega_p1.'*(I_p1*omega_p1));
T_p2 = simplify(1/2*m_p*VG2.'*VG2 + 1/2*omega_p2.'*(I_p2*omega_p2));
T = simplify(T_s + T_p1 + T_p2); 
L = T; %Since the potential energy is zero, L = T
%% Equations of motion
q = [q1;q2;q3;q4;q5;q6;q7;q8];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7;dq8]; 
Q = jacobian(omega_s, dq).'*[tau1; tau2; tau3] + jacobian(omega_p1, dq).'*[0;0;tau_p1] + jacobian(omega_p2, dq).'*[0;0;tau_p2];
Q = simplify(Q);

pL_pq = jacobian(L,q).';
pL_pdq = jacobian(L,dq).';
M = jacobian(pL_pdq,dq);
M = simplify(M);
B = jacobian(pL_pdq,q)*dq + diff(pL_pdq, t) - pL_pq;
B = simplify(B);

%% Linear momentum
Ptot = simplify(m_s*VG_XYZ + m_p*(VG1_XYZ + VG2_XYZ));

%% Functions production
matlabFunction(T,'File','Energy');
matlabFunction(Ptot,'File','Lmoment');
matlabFunction(M,'File','MassMat');
matlabFunction(B,'File','Bias');
matlabFunction(Q,'File','GenForc');