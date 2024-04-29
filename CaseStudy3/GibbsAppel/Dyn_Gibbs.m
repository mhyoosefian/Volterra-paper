function dz = Dyn_Gibbs(t,z) %#ok<*INUSL>
dz = z;
global m_sat m Ixx Iyy Izz Ixy Ixz Iyz l a tau1 tau2 tau3
F = -0.018*sin(0.089*t) + 0.012*cos(0.0485*t);
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); q7 = z(7);
U1 = z(8); U2 = z(9); U3 = z(10); U4 = z(11); U5 = z(12); U6 = z(13); U7 = z(14);
dq = dq_itu_Gibbs(U1,U2,U3,U4,U5,U6,U7,q2,q3);
M_g = M_gibss(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,l,m,m_sat,q1,q2,q3,q4);
B_g = B_gibss(F,Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,a,l,m,q1,q2,q3,q4,tau1,tau2,tau3);
POW = Power(F,U1,U2,U3,U4,q2,q3,tau1,tau2,tau3);
dz = [dq; -inv(M_g)*B_g; POW];