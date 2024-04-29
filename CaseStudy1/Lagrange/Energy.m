function E = Energy(dq1,dq2,dq3,g,l,m1,m2,q1,q2)
%ENERGY
%    E = ENERGY(DQ1,DQ2,DQ3,G,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    02-Oct-2020 17:51:13

t2 = dq3.^2;
t3 = l.^2;
t4 = cos(q1);
t5 = cos(q2);
E = m1.*t2.*(1.0./2.0)+m2.*t2+dq1.^2.*m2.*t3.*(2.0./3.0)+dq2.^2.*m2.*t3.*(1.0./6.0)-g.*l.*m1-g.*l.*m2.*t4.*(3.0./2.0)-g.*l.*m2.*t5.*(1.0./2.0)+dq1.*dq2.*m2.*t3.*cos(q1-q2).*(1.0./2.0)+dq1.*dq3.*l.*m2.*t4.*(3.0./2.0)+dq2.*dq3.*l.*m2.*t5.*(1.0./2.0);