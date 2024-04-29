function dq_itu = dq_itu_Voltra(U1,U2,U3,U4,U5,b,m_p,m_s,px,py,pz,q1,q2,q3,q4,q5)
%DQ_ITU_VOLTRA
%    DQ_ITU = DQ_ITU_VOLTRA(U1,U2,U3,U4,U5,B,M_P,M_S,PX,PY,PZ,Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    26-Aug-2020 09:59:23

t2 = cos(q3);
t3 = sin(q3);
t4 = cos(q2);
t5 = 1.0./t4;
t6 = sin(q2);
t7 = m_p.*2.0;
t8 = m_s+t7;
t9 = 1.0./t8;
t10 = cos(q1);
t11 = cos(q4);
t12 = sin(q1);
t13 = cos(q5);
t14 = sin(q4);
t15 = t4.*t10.*t14;
t16 = sin(q5);
t17 = t4.*t10.*t16;
t18 = t2.*t11.*t12;
t19 = t3.*t6.*t10.*t13;
t20 = m_p.*4.0;
t21 = m_s.*2.0;
t22 = t20+t21;
t23 = 1.0./t22;
t24 = t3.*t12;
t25 = t2.*t6.*t10;
t26 = t24+t25;
t27 = t2.*t10.*t13;
t28 = t4.*t12.*t16;
t29 = t3.*t6.*t12.*t13;
t30 = t2.*t10.*t11;
t31 = t3.*t6.*t11.*t12;
t32 = t14-t16;
t33 = t3.*t10;
t34 = t33-t2.*t6.*t12;
t35 = t11+t13-2.0;
t36 = t6.*t14;
t37 = t3.*t4.*t11;
t38 = t3.*t4.*t13;
dq_itu = [t5.*(U2.*t3+U3.*t2);U2.*t2-U3.*t3;t5.*(U1.*t4+U2.*t3.*t6+U3.*t2.*t6);U4;U5;px.*t9-U4.*b.*m_p.*t9.*(t15+t18-t3.*t6.*t10.*t11).*(1.0./2.0)-U3.*b.*m_p.*t23.*(t15-t17+t18-t19-t2.*t12.*2.0+t3.*t6.*t10.*2.0+t2.*t12.*t13-t3.*t6.*t10.*t11)-U5.*b.*m_p.*t9.*(t17+t19-t2.*t12.*t13).*(1.0./2.0)-U2.*b.*m_p.*t9.*t26.*t35.*(1.0./2.0)+U1.*b.*m_p.*t23.*t26.*t32;py.*t9+U3.*b.*m_p.*t23.*(t27+t28+t29+t30+t31-t2.*t10.*2.0-t3.*t6.*t12.*2.0-t4.*t12.*t14)-U5.*b.*m_p.*t9.*(t27+t28+t29).*(1.0./2.0)+U4.*b.*m_p.*t23.*(t30+t31-t4.*t12.*t14)-U1.*b.*m_p.*t9.*t32.*t34.*(1.0./2.0)+U2.*b.*m_p.*t9.*t34.*t35.*(1.0./2.0);pz.*t9-U5.*b.*m_p.*t9.*(t38-t6.*t16).*(1.0./2.0)+U4.*b.*m_p.*t9.*(t36+t37).*(1.0./2.0)+U3.*b.*m_p.*t23.*(t36+t37+t38-t3.*t4.*2.0-t6.*t16)-U2.*b.*m_p.*t2.*t4.*t9.*t35.*(1.0./2.0)+U1.*b.*m_p.*t2.*t4.*t23.*t32];