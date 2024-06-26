function dq_itu = dq_itu_Voltra(P,U,l,m1,m2,q1,q2)
%DQ_ITU_VOLTRA
%    DQ_ITU = DQ_ITU_VOLTRA(P,U,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    03-Oct-2020 21:52:28

t2 = q1-q2;
t3 = cos(t2);
t4 = t3+1.0;
t5 = 1.0./t4;
t6 = m2.*2.0;
t7 = m1+t6;
t8 = 1.0./t7;
dq_itu = [-U.*t5;U.*t3.*t5;P.*t8-U.*l.*m2.*t5.*t8.*(cos(q1-q2.*2.0)-cos(q1).*5.0).*(1.0./4.0)];
