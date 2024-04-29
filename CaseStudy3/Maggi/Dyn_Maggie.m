function dz = Dyn_Maggie(t,z) 
dz = z;
global m_sat m Ixx Iyy Izz Ixy Ixz Iyz l a tau1 tau2 tau3
F = -0.018*sin(0.089*t) + 0.012*cos(0.0485*t);
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); q7 = z(7);
dq1 = z(8); dq2 = z(9); dq3 = z(10); dq4 = z(11); dq5 = z(12); dq6 = z(13); dq7 = z(14);
dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
B_ = B(F,Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,dq1,dq2,dq3,dq4,l,m,q1,q2,q3,q4,tau1,tau2,tau3);
A_ = A(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,l,m,m_sat,q1,q2,q3,q4);
POW = Power(F,dq1,dq2,dq3,dq4,q2,q3,tau1,tau2,tau3);
dz = [dq; -inv(A_)*B_; POW];
