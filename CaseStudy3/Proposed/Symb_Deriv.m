clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi   s   X    Y    Z  
syms q1     q2      q3    q4  q5   q6   q7  real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real
syms m_sat m Ixx Iyy Izz Ixy Ixz Iyz l a real
syms tau1 tau2 tau3 F real
syms U1 U2 U3 U4 real
syms dU1 dU2 dU3 dU4 real
syms px py pz real

P = [px; py; pz];

U = [U1; U2; U3; U4];
dU = [dU1; dU2; dU3; dU4];
q = [q1;q2;q3;q4;q5;q6;q7];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7]; 
dq_I = [dq5;dq6;dq7];
dq_NI = [dq1;dq2;dq3;dq4];

Rphi = [1 0 0;0 cos(q3) sin(q3);0 -sin(q3) cos(q3)];
Rtheta = [cos(q2) 0 -sin(q2);0 1 0;sin(q2) 0 cos(q2)];
Rpsi = [cos(q1) sin(q1) 0;-sin(q1) cos(q1) 0;0 0 1];
R = simplify(Rphi*Rtheta*Rpsi);

%% Angular velocities
omega_s = simplify(Rphi*Rtheta*[0;0;dq1] + Rphi*[0;dq2;0] + [dq3;0;0]); %expressed in xyz(satelite coordinate system)

%% Linear velocities
VG_XYZ = [dq5;dq6;dq7];
VG = simplify(R*VG_XYZ); %expressed in xyz(satelite coordinate system)
Vm = simplify(VG + cross(omega_s, [0;0;l/2+a+q4]) + [0;0;dq4]);
Vm_XYZ = simplify(R.'*Vm);

PXYZ = simplify(m_sat*VG_XYZ + m*Vm_XYZ);

%% Moments of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)

%% Kinetic Energy
T_s = simplify(1/2*m_sat*VG.'*VG+1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_m = simplify(1/2*m*Vm.'*Vm);
T = simplify(T_s+T_m);

%% Quasi-velocities
pT_pdq = jacobian(T, dq).';
M = (jacobian(pT_pdq, dq));
M11 = M(1:4,1:4);
M12 = M(1:4,5:7);
M21 = M12.';
M22 = M(5:7,5:7);
M2 = simplify([M21 M22]);
N2 = simplify(jacobian(T, dq_I).' - M21*dq_NI - M22*dq_I);

u1 = omega_s(1);
u2 = omega_s(2);
u3 = omega_s(3);
% u1 = dq1;
% u2 = dq2;
% u3 = dq3;
u4 = dq4;
u = [u1; u2; u3; u4];
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = simplify(inv([Y; M2])*[eye(4,4); zeros(3,4)]);
X = simplify(-inv([Y; M2])*[Z; N2-P]);
dq_itu = simplify(W*U+X);

%% Generalized Forces
for i=1:length(q)
Q(i,:) = dot(jacobian(omega_s, dq(i)), [tau1; tau2; tau3]) + dot(jacobian(VG, dq(i)), [0; 0; -F])...
    + dot(jacobian(Vm, dq(i)), [0; 0; F]);
end
Qpr = W.'*Q;

%% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VG = simplify(subs(VG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
Vm = simplify(subs(Vm, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
VG_XYZ = simplify(subs(VG_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
Vm_XYZ = simplify(subs(Vm_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
omega_s = simplify(subs(omega_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

%% Kinetic Energy 
T = simplify(subs(T, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

%% Reduced Equations of Voltra
pT_pU = simplify(jacobian(T, U).');
d_dt_pT_puNI = zeros(size(pT_pU));
for i =1:length(q)
    d_dt_pT_puNI = d_dt_pT_puNI + diff(pT_pU, q(i))*dq_itu(i);
end
d_dt_pT_puNI = (d_dt_pT_puNI);
for i =1:length(U)
    d_dt_pT_puNI = d_dt_pT_puNI + diff(pT_pU, U(i))*dU(i);
end
d_dt_pT_puNI = (subs(d_dt_pT_puNI, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

pVG_pU = jacobian(VG, U);
pVm_pU = jacobian(Vm, U);
d_dt_pVG_pU = zeros(size(pVG_pU));
d_dt_pVm_pU = zeros(size(pVm_pU));
for i =1:length(q)
    d_dt_pVG_pU = d_dt_pVG_pU + diff(pVG_pU, q(i))*dq_itu(i);
    d_dt_pVm_pU = d_dt_pVm_pU + diff(pVm_pU, q(i))*dq_itu(i);
end
d_dt_pVG_pU = simplify(d_dt_pVG_pU);
d_dt_pVm_pU = simplify(d_dt_pVm_pU);
for i =1:length(U)
    d_dt_pVG_pU = d_dt_pVG_pU + diff(pVG_pU, U(i))*dU(i);
    d_dt_pVm_pU = d_dt_pVm_pU + diff(pVm_pU, U(i))*dU(i);
end
d_dt_pVG_pU = (subs(d_dt_pVG_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
d_dt_pVm_pU = (subs(d_dt_pVm_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

for i=1:length(U)
    d_dt_pVG_pU(:,i) = (d_dt_pVG_pU(:,i) + cross(omega_s, pVG_pU(:,i)));
    d_dt_pVm_pU(:,i) = (d_dt_pVm_pU(:,i) + cross(omega_s, pVm_pU(:,i)));
end
PG = m_sat*VG;
Pm = m*Vm;


pomegas_pU = jacobian(omega_s, U);
d_dt_pomegas_pU = zeros(size(pomegas_pU));
for i =1:length(q)
    d_dt_pomegas_pU = d_dt_pomegas_pU + diff(pomegas_pU, q(i))*dq_itu(i);
end
d_dt_pomegas_pU = simplify(d_dt_pomegas_pU);

for i =1:length(U)
    d_dt_pomegas_pU = d_dt_pomegas_pU + diff(pomegas_pU, U(i))*dU(i);
end
d_dt_pomegas_pU = (subs(d_dt_pomegas_pU, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

for i=1:length(U)
    d_dt_pomegas_pU(:,i) = (d_dt_pomegas_pU(:,i) + cross(omega_s, pomegas_pU(:,i)));
end

HG = I_s*omega_s; %expressed in xyz(satelite coordinate system)


eq = d_dt_pT_puNI - simplify(d_dt_pVG_pU.'*PG) - (d_dt_pVm_pU.'*Pm)...
- (d_dt_pomegas_pU.'*HG) - Qpr;
eq = simplify(eq);
A_ = simplify(jacobian(eq, dU));
B_ = simplify(eq - A_*dU);

%% Linear momentum
Ptot = (m_sat*VG_XYZ+m*Vm_XYZ);

%% Power
% Pow = dot([0; 0; -F], VG) + dot([0; 0; F], Vm) + dot([tau1; tau2; tau3], omega_s);
Pow = Qpr.'*U;
%% Functions production
matlabFunction(T,'File','Energy');
matlabFunction(Ptot,'File','Lmoment');
matlabFunction(A_,'File','MassMat');
matlabFunction(B_,'File','Bias');
matlabFunction(dq_itu,'File','dq_itu_Voltra');
matlabFunction([u1; u2; u3; u4],'File','U_NI');     %Is used for initial state vector
matlabFunction(PXYZ,'File','PXYZ_Voltra');
matlabFunction(M21,'File','M21_Voltra');
matlabFunction(M22,'File','M22_Voltra');
matlabFunction(Pow,'File','Power');