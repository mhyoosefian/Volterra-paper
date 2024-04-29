function dz = Dyn_Maggie(t,z) %#ok<*INUSL>
global m_s m_p a b c Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 tau_p1 tau_p2
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); q7 = z(7); q8 = z(8); %#ok<*NASGU>
dq1 = z(9); dq2 = z(10); dq3 = z(11); dq4 = z(12); dq5 = z(13); dq6 = z(14); dq7 = z(15); dq8 = z(16); %#ok<*NASGU>
dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7; dq8];
B_ = B(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,b,c,dq1,dq2,dq3,dq4,dq5,m_p,q1,q2,q3,q4,q5,tau1,tau2,tau3,tau_p1,tau_p2);
A_ = A(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,b,c,m_p,m_s,q1,q2,q3,q4,q5);
ddq = -inv(A_)*B_;
dz = [dq; ddq];