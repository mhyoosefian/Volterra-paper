function dz = Dyn_Gibbs(t,z) %#ok<*INUSL>
% dz = z;
global m_s m_p a b c Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 tau_p1 tau_p2

q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5);
U1 = z(9); U2 = z(10); U3 = z(11); U4 = z(12); U5 = z(13); U6 = z(14); U7 = z(15); U8 = z(16);
dq = dq_itu_Gibbs(U1,U2,U3,U4,U5,U6,U7,U8,q2,q3);
M_g = M_gibss(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,b,c,m_p,m_s,q1,q2,q3,q4,q5);
B_g = B_gibss(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,U5,a,b,c,m_p,q1,q2,q3,q4,q5,tau1,tau2,tau3,tau_p1,tau_p2);
% du = -inv(M_g)*B_g;
dz = [dq; -inv(M_g)*B_g];