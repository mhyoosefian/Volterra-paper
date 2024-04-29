function dz = Dyn_LagIM(t,z) %#ok<*INUSL>

global m1 l m2 g tau

q1 = z(1); q2 = z(2); q3 = z(3); %#ok<*NASGU>
dq1 = z(4); dq2 = z(5); dq3 = z(6);

B_IM = B(dq1,dq2,dq3,g,l,m2,q1,q2,tau);
A_IM = A(l,m1,m2,q1,q2);
dz = inv(A_IM)*B_IM;