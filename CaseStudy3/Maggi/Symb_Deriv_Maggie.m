clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi   s   X    Y    Z  
syms q1     q2      q3    q4  q5   q6   q7  real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real
syms m_sat m Ixx Iyy Izz Ixy Ixz Iyz l a real
syms tau1 tau2 tau3 F real

syms U1 U2 U3 U4 U5 U6 U7 real

U = [U1; U2; U3; U4; U5; U6; U7];
q = [q1;q2;q3;q4;q5;q6;q7];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7]; 
ddq = [ddq1;ddq2;ddq3;ddq4;ddq5;ddq6;ddq7]; 

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

%% Moments of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)

%% Kinetic energy
T_s = simplify(1/2*m_sat*VG.'*VG+1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_m = simplify(1/2*m*Vm.'*Vm);
T = simplify(T_s+T_m); 

%% W Matrix and Quasi-velocities
u1 = omega_s(1);
u2 = omega_s(2);
u3 = omega_s(3);
u4 = dq4;
u5 = dq5;
u6 = dq6;
u7 = dq7;
u = [u1; u2; u3; u4; u5; u6; u7];
Y = jacobian(u,dq);
Z = simplify(u - Y*dq);
W = simplify(inv(Y));
X = simplify(-Y\Z);

%% Generalized Forces
for i=1:length(q)
Q(i,:) = dot(jacobian(omega_s, dq(i)), [tau1; tau2; tau3]) + dot(jacobian(VG, dq(i)), [0; 0; -F])...
    + dot(jacobian(Vm, dq(i)), [0; 0; F]);
end
Qpr = W.'*Q;
%% Maggies' Equations
pT_pq = jacobian(T, q).';
pT_pdq = jacobian(T, dq).';
d_dt_pT_pdq = zeros(size(dq));
for i=1:length(q)
   d_dt_pT_pdq =  d_dt_pT_pdq + diff(pT_pdq, dq(i))*ddq(i) + + diff(pT_pdq, q(i))*dq(i);
end
eq = simplify(W.'*(d_dt_pT_pdq - pT_pq) - Qpr);
M_ = jacobian(eq, ddq);
B_ = simplify(eq - M_*ddq);
Ptot = (m_sat*VG_XYZ+m*Vm_XYZ);
%% Power
Pow = Q.'*dq;

%% Functions
matlabFunction(M_, 'File', 'A');
matlabFunction(B_,'File','B');
matlabFunction(T, 'File', 'Energy_Maggie');
matlabFunction(Ptot, 'File', 'Lmom_Maggie');
matlabFunction(Pow,'File','Power');