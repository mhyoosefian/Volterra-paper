function M_g = M_gibss(l,m1,m2,q1,q2)
%M_GIBSS
%    M_G = M_GIBSS(L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    04-Oct-2020 20:58:07

t2 = q1-q2;
t3 = cos(t2);
t4 = q1.*2.0;
t19 = q2.*2.0;
t5 = t4-t19;
t6 = cos(t5);
t7 = q1.*3.0;
t8 = q2.*-3.0+t7;
t9 = cos(t8);
t10 = q1.*4.0;
t11 = q2.*-4.0+t10;
t12 = cos(t11);
t13 = q1.*5.0;
t14 = q2.*-5.0+t13;
t15 = cos(t14);
t16 = q1.*6.0;
t17 = q2.*-6.0+t16;
t18 = cos(t17);
t20 = t3+1.0;
t21 = 1.0./t20;
t22 = q1-t19;
t23 = cos(t22);
t24 = cos(q1);
t25 = t23-t24.*5.0;
M_g = reshape([(l.^2.*m2.*(t3.*-2.72e2-t6.*9.7e1+t9.*8.0+t12.*2.2e1+t15.*8.0+t18-1.82e2).*(-2.0./3.0))./(t3.*7.92e2+t6.*4.95e2+t9.*2.2e2+t12.*6.6e1+t15.*1.2e1+t18+4.62e2),l.*m2.*t21.*t25.*(-1.0./4.0),l.*m2.*t21.*t25.*(-1.0./4.0),m1+m2.*2.0],[2,2]);