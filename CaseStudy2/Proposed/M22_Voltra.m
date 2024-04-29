function M22 = M22_Voltra(m_p,m_s,q1,q2,q3)
%M22_VOLTRA
%    M22 = M22_VOLTRA(M_P,M_S,Q1,Q2,Q3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    26-Aug-2020 09:59:26

t3 = cos(q3);
t4 = sin(q1);
t5 = cos(q1);
t6 = sin(q2);
t7 = sin(q3);
t9 = t4.*t7;
t10 = t3.*t5.*t6;
t2 = t9+t10;
t12 = t3.*t4;
t13 = t5.*t6.*t7;
t8 = t12-t13;
t11 = t2.^2;
t14 = t8.^2;
t15 = cos(q2);
t16 = t5.^2;
t17 = t15.^2;
t18 = t3.*t5;
t19 = t4.*t6.*t7;
t20 = t18+t19;
t21 = t5.*t7;
t23 = t3.*t4.*t6;
t22 = t21-t23;
t24 = m_p.*t4.*t5.*t17.*2.0;
t25 = m_s.*t4.*t5.*t17;
t26 = t24+t25-m_p.*t2.*t22.*2.0-m_p.*t8.*t20.*2.0-m_s.*t2.*t22-m_s.*t8.*t20;
t27 = t20.^2;
t28 = t22.^2;
t29 = t4.^2;
t30 = m_p.*t2.*t3.*t15.*2.0;
t31 = m_s.*t2.*t3.*t15;
t32 = t30+t31-m_p.*t5.*t6.*t15.*2.0-m_p.*t7.*t8.*t15.*2.0-m_s.*t5.*t6.*t15-m_s.*t7.*t8.*t15;
t33 = m_p.*t7.*t15.*t20.*2.0;
t34 = m_s.*t7.*t15.*t20;
t35 = t33+t34-m_p.*t4.*t6.*t15.*2.0-m_p.*t3.*t15.*t22.*2.0-m_s.*t4.*t6.*t15-m_s.*t3.*t15.*t22;
t36 = t6.^2;
t37 = t3.^2;
t38 = t7.^2;
M22 = reshape([m_p.*t11.*2.0+m_p.*t14.*2.0+m_s.*t11+m_s.*t14+m_p.*t16.*t17.*2.0+m_s.*t16.*t17,t26,t32,t26,m_p.*t27.*2.0+m_p.*t28.*2.0+m_s.*t27+m_s.*t28+m_p.*t17.*t29.*2.0+m_s.*t17.*t29,t35,t32,t35,m_p.*t36.*2.0+m_s.*t36+m_p.*t17.*t37.*2.0+m_p.*t17.*t38.*2.0+m_s.*t17.*t37+m_s.*t17.*t38],[3,3]);
