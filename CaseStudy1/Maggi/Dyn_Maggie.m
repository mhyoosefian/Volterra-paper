function dz = Dyn_Maggie(t,z) %#ok<*INUSL>

global m1 l m2 g tau

q1 = z(1); q2 = z(2); q3 = z(3); %#ok<*NASGU>
dq1 = z(4); dq2 = z(5); dq3 = z(6);
dq = [dq1; dq2; dq3];
% ddq1 = 0; ddq2 = 0; ddq3 = 0;
B_ = B_Maggie(dq1,dq2,g,l,m2,q1,q2,tau);
A_ = A_Maggie(l,m1,m2,q1,q2);
% ddq = -inv(A_)*B_;
dz = [dq; -inv(A_)*B_];