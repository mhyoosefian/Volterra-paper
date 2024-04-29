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
syms dU1 dU2 dU3 dU4 dU5 dU6 dU7 real

U = [U1; U2; U3; U4; U5; U6; U7];
dU = [dU1; dU2; dU3; dU4; dU5; dU6; dU7];
q = [q1; q2; q3; q4; q5; q6; q7];
dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7]; 

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

%% Quasi-velocities
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
dq_itu = W*U+X;

%% Generalized Forces
for i=1:length(q)
Q(i,:) = dot(jacobian(omega_s, dq(i)), [tau1; tau2; tau3]) + dot(jacobian(VG, dq(i)), [0; 0; -F])...
    + dot(jacobian(Vm, dq(i)), [0; 0; F]);
end
Qpr = W.'*Q;

%% Stating Velocites in terms of u                  %Velocities presented here are funcitons of q and U
VG = (subs(VG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
VG_XYZ = (subs(VG_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
Vm = (subs(Vm, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
Vm_XYZ = (subs(Vm_XYZ, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));
omega_s = (subs(omega_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

%% Kinetic energy
T_s = simplify(1/2*m_sat*VG.'*VG+1/2*omega_s.'*(I_s*omega_s)); %#ok<*MHERM>
T_m = simplify(1/2*m*Vm.'*Vm);
T = simplify(T_s+T_m);

%% Linear accelerations
accelG = zeros(3,1);
for i=1:length(q)
accelG = accelG + diff(VG,q(i))*dq(i);
end
for i=1:length(U)
accelG = accelG + diff(VG,U(i))*dU(i);
end
accelG = simplify(accelG + cross(omega_s, VG));
accelG = (subs(accelG, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

accelm = zeros(3,1);
for i=1:length(q)
accelm = accelm + diff(Vm,q(i))*dq(i);
end
for i=1:length(U)
accelm = accelm + diff(Vm,U(i))*dU(i);
end
accelm = simplify(accelm + cross(omega_s, Vm));
accelm = (subs(accelm, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

%% Angular accelerations
alpha_s = zeros(3,1);
for i=1:length(q)
alpha_s = alpha_s + diff(omega_s,q(i))*dq(i);
end
for i=1:length(U)
alpha_s = alpha_s + diff(omega_s,U(i))*dU(i);
end
alpha_s = simplify(alpha_s);
alpha_s = (subs(alpha_s, {dq1 dq2 dq3 dq4 dq5 dq6 dq7}, {dq_itu(1) dq_itu(2) dq_itu(3) dq_itu(4) dq_itu(5) dq_itu(6) dq_itu(7)}));

%% Building Gibss Function
S = 0.5*(m_sat*(dot(accelG,accelG)) + m*(dot(accelm,accelm)))...
    + 0.5*(dot(alpha_s,I_s*alpha_s))...
    + dot(alpha_s,cross(omega_s, I_s*omega_s));

%% Equations of Motion
Sdot1 = diff(S,dU1);
Sdot2 = diff(S,dU2);
Sdot3 = diff(S,dU3);
Sdot4 = diff(S,dU4);
Sdot5 = diff(S,dU5);
Sdot6 = diff(S,dU6);
Sdot7 = diff(S,dU7);

eq = [Sdot1; Sdot2; Sdot3; Sdot4; Sdot5; Sdot6; Sdot7] - Qpr;
M_g = simplify(jacobian(eq, dU));
B_g = simplify(eq - M_g*dU);

Ptot = (m_sat*VG_XYZ+m*Vm_XYZ);

%% Power
Pow = Qpr.'*U;

%% Functions
matlabFunction(M_g, 'File', 'M_gibss');
matlabFunction(B_g, 'File', 'B_gibss');
matlabFunction(T, 'File', 'Energy_Gibbs');
matlabFunction(Ptot, 'File', 'ptot_Gibbs');
matlabFunction(dq_itu, 'File', 'dq_itu_Gibbs');
matlabFunction(u, 'File', 'u_Gibbs');     % Will be used in initial condition
matlabFunction(Pow,'File','Power');