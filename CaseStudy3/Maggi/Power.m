function Pow = Power(F,dq1,dq2,dq3,dq4,q2,q3,tau1,tau2,tau3)
%POWER
%    POW = POWER(F,DQ1,DQ2,DQ3,DQ4,Q2,Q3,TAU1,TAU2,TAU3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    07-Oct-2020 06:47:14

t2 = cos(q3);
t3 = cos(q2);
t4 = sin(q3);
Pow = F.*dq4+dq3.*tau1+dq1.*(-tau1.*sin(q2)+t2.*t3.*tau3+t3.*t4.*tau2)+dq2.*(t2.*tau2-t4.*tau3);
