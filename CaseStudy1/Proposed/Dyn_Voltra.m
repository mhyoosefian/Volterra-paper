function dz = Dyn_Voltra(t,z)
global m1 m2 l tau g P
q1 = z(1); q2 = z(2);
U = z(4);
B_ = B_Voltra(U,g,l,m1,m2,q1,q2,tau);
A_ = A_Voltra(l,m1,m2,q1,q2);
dq = dq_itu_Voltra(P,U,l,m1,m2,q1,q2);
dz = [dq; -inv(A_)*B_];