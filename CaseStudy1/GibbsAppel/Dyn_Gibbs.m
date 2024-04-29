function dz = Dyn_Gibbs(t,z) %#ok<*INUSL>
global m1 l m2 g tau
q1 = z(1); q2 = z(2);
U1 = z(4); U2 = z(5);
dU1 = 0; dU2 = 0;
% dq = [U1/(cos(q1 - q2) + 1); -(U1*cos(q1 - q2))/(cos(q1 - q2) + 1); U2];
dq = dq_itu_Gibbs(U1,U2,q1,q2);
M_g = M_gibss(l,m1,m2,q1,q2);
B_g = B_gibss(U1,U2,dU1,dU2,g,l,m1,m2,q1,q2,tau);
dz = [dq; -inv(M_g)*B_g];