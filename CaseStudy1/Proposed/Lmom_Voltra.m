function Ptot = Lmom_Voltra(P,U,l,m1,m2,q1,q2)
%LMOM_VOLTRA
%    PTOT = LMOM_VOLTRA(P,U,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    03-Oct-2020 21:52:28

t2 = q1-q2;
t3 = cos(t2);
t4 = t3+1.0;
t5 = 1.0./t4;
t6 = cos(q1);
t7 = m2.*2.0;
t8 = m1+t7;
t9 = 1.0./t8;
t10 = q1-q2.*2.0;
t11 = cos(t10);
t12 = t6.*5.0;
t13 = t11-t12;
Ptot = [m2.*(P.*t9.*2.0-U.*l.*t5.*t6.*(3.0./2.0)+U.*l.*t3.*t5.*cos(q2).*(1.0./2.0)-U.*l.*m2.*t5.*t9.*t13.*(1.0./2.0))+m1.*(P.*t9-U.*l.*m2.*t5.*t9.*t13.*(1.0./4.0));-m2.*(U.*l.*t5.*sin(q1).*(3.0./2.0)-U.*l.*t3.*t5.*sin(q2).*(1.0./2.0));0.0];
