function B_ = B_Maggie(dq1,dq2,g,l,m2,q1,q2,tau)
%B_MAGGIE
%    B_ = B_MAGGIE(DQ1,DQ2,G,L,M2,Q1,Q2,TAU)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    05-Oct-2020 08:20:54

t2 = q1-q2;
t3 = cos(t2);
t4 = l.^2;
t5 = dq1.^2;
t6 = sin(q1);
t7 = dq2.^2;
t8 = sin(t2);
B_ = [(-tau-t3.*tau+g.*l.*m2.*t6.*(5.0./4.0)+m2.*t4.*t5.*sin(q1.*2.0-q2.*2.0).*(1.0./4.0)+m2.*t4.*t7.*t8.*(1.0./2.0)+g.*l.*m2.*sin(q1-q2.*2.0).*(1.0./4.0))./(t3+1.0);l.*m2.*(t5.*t6.*3.0+t7.*sin(q2)).*(-1.0./2.0);-dq1.*(dq1.*l.*t8-dq2.*l.*t8)];
