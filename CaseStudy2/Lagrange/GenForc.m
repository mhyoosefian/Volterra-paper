function Q = GenForc(q2,q3,tau1,tau2,tau3,tau_p1,tau_p2)
%GENFORC
%    Q = GENFORC(Q2,Q3,TAU1,TAU2,TAU3,TAU_P1,TAU_P2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    21-Aug-2020 12:35:02

t2 = cos(q2);
t3 = cos(q3);
t4 = sin(q3);
Q = [-tau1.*sin(q2)+t2.*t3.*tau3+t2.*t4.*tau2+t2.*t3.*tau_p1+t2.*t3.*tau_p2;t3.*tau2-t4.*tau3-t4.*tau_p1-t4.*tau_p2;tau1;tau_p1;-tau_p2;0.0;0.0;0.0];