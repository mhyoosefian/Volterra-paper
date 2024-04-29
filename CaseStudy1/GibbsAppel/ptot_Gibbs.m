function Ptot = ptot_Gibbs(dq1,dq2,dq3,l,m1,m2,q1,q2)
%PTOT_GIBBS
%    PTOT = PTOT_GIBBS(DQ1,DQ2,DQ3,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    04-Oct-2020 20:58:09

t2 = cos(q1);
t3 = sin(q1);
Ptot = [dq3.*m1+m2.*(dq3+dq1.*l.*t2.*(1.0./2.0))+m2.*(dq3+dq2.*l.*cos(q2).*(1.0./2.0)+dq1.*l.*t2);m2.*(dq2.*l.*sin(q2).*(1.0./2.0)+dq1.*l.*t3)+dq1.*l.*m2.*t3.*(1.0./2.0);0.0];
