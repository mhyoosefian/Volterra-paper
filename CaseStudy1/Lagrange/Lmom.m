function Ptot = Lmom(dq1,dq2,dq3,l,m1,m2,q1,q2)
%LMOM
%    PTOT = LMOM(DQ1,DQ2,DQ3,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    02-Oct-2020 17:51:13

Ptot = [dq3.*m1+m2.*(dq3.*2.0+dq1.*l.*cos(q1).*(3.0./2.0)+dq2.*l.*cos(q2).*(1.0./2.0));m2.*(dq1.*l.*sin(q1).*(3.0./2.0)+dq2.*l.*sin(q2).*(1.0./2.0));0.0];
