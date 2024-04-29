clc;
clear all %#ok<*CLALL>

syms t
     %psi  theta    phi   xi   X    Y    Z  
syms q1     q2      q3    q4  q5   q6   q7  real 
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real
syms m_sat m Ixx Iyy Izz Ixy Ixz Iyz l a real
syms tau1 tau2 tau3 F real

Rphi = [1 0 0;0 cos(q3) sin(q3);0 -sin(q3) cos(q3)];
Rtheta = [cos(q2) 0 -sin(q2);0 1 0;sin(q2) 0 cos(q2)];
Rpsi = [cos(q1) sin(q1) 0;-sin(q1) cos(q1) 0;0 0 1];
R = simplify(Rphi*Rtheta*Rpsi);


%% Angular velocities
omega_s = simplify(Rphi*Rtheta*[0;0;dq1] + Rphi*[0;dq2;0] + [dq3;0;0]); %expressed in xyz(satelite coordinate system)

%% Velocities of centers of mass
VG_XYZ = [dq5;dq6;dq7];
VG = simplify(R*VG_XYZ); %expressed in xyz(satelite coordinate system)

Vm = simplify(VG + cross(omega_s, [0;0;l/2+a+q4]) + [0;0;dq4]);
Vm_XYZ = simplify(R.'*Vm);

%% Moment of Inertia
I_s = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz]; %expressed in xyz(satelite coordinate system)

%% Lagrangian function
T_s = simplify(1/2*m_sat*VG.'*VG+1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_m = simplify(1/2*m*Vm.'*Vm);
T = simplify(T_s+T_m); 

%% Equations of motion
q = [q1;q2;q3;q4;q5;q6;q7];
dq = [dq1;dq2;dq3;dq4;dq5;dq6;dq7]; 
ddq = [ddq1;ddq2;ddq3;ddq4;ddq5;ddq6;ddq7]; 
L = T; %Since the potential energy is zero, L = T

for i=1:length(q)
Q(i,:) = dot(jacobian(omega_s, dq(i)), [tau1; tau2; tau3]) + dot(jacobian(VG, dq(i)), [0; 0; -F])...
    + dot(jacobian(Vm, dq(i)), [0; 0; F]);
end

pL_pq = jacobian(L,q).';
pL_pdq = jacobian(L,dq).';
M = jacobian(pL_pdq,dq);
B = jacobian(pL_pdq,q)*dq + diff(pL_pdq, t) - pL_pq;

%% Linear momentum
Ptot = (m_sat*VG_XYZ+m*Vm_XYZ);

%% Power
Pow = Q.'*dq;
%% Functions production
matlabFunction(B,'File','B');
matlabFunction(M,'File','M');
matlabFunction(T,'File','Energy');
matlabFunction(Ptot,'File','Lmoment');
matlabFunction(Q,'File','GenForce');
matlabFunction(Pow,'File','Power');