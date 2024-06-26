function u = u_Gibbs(dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,q2,q3)
%U_GIBBS
%    U = U_GIBBS(DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,DQ8,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    30-Sep-2020 20:13:32

t2 = sin(q3);
t3 = cos(q2);
t4 = cos(q3);
u = [dq3-dq1.*sin(q2);dq2.*t4+dq1.*t2.*t3;-dq2.*t2+dq1.*t3.*t4;dq4;dq5;dq6;dq7;dq8];
