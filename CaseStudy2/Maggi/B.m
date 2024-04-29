function B_ = B(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,b,c,dq1,dq2,dq3,dq4,dq5,m_p,q1,q2,q3,q4,q5,tau1,tau2,tau3,tau_p1,tau_p2)
%B
%    B_ = B(IXX,IXY,IXZ,IYY,IYZ,IZZ,A,B,C,DQ1,DQ2,DQ3,DQ4,DQ5,M_P,Q1,Q2,Q3,Q4,Q5,TAU1,TAU2,TAU3,TAU_P1,TAU_P2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    23-Aug-2020 11:53:15

t2 = cos(q2);
t3 = dq2.^2;
t4 = cos(q3);
t5 = q3.*2.0;
t6 = sin(t5);
t7 = dq1.^2;
t8 = t2.^2;
t9 = t4.^2;
t10 = a.^2;
t11 = sin(q2);
t12 = b.^2;
t13 = sin(q3);
t14 = cos(q4);
t15 = cos(q5);
t16 = c.^2;
t17 = t14.^2;
t18 = t15.^2;
t19 = sin(q4);
t20 = sin(q5);
t21 = dq3.^2;
t22 = q4.*2.0;
t23 = sin(t22);
t24 = q5.*2.0;
t25 = sin(t24);
t26 = dq4.^2;
t27 = dq5.^2;
t28 = m_p.*t7.*t8.*t12.*t14.*t19.*(2.0./3.0);
t29 = m_p.*t3.*t9.*t12.*t14.*t19.*(1.0./3.0);
t30 = a.*b.*m_p.*t3.*t9.*t14.*(1.0./4.0);
t31 = dq2.*dq3.*m_p.*t4.*t12.*t17.*(2.0./3.0);
t32 = a.*b.*dq1.*dq3.*m_p.*t11.*t14.*(1.0./2.0);
t33 = dq1.*dq3.*m_p.*t2.*t12.*t13.*t17.*(2.0./3.0);
t34 = dq1.*dq3.*m_p.*t11.*t12.*t14.*t19.*(2.0./3.0);
t35 = dq1.*dq2.*m_p.*t2.*t4.*t12.*t13.*t14.*t19.*(2.0./3.0);
t36 = a.*b.*dq1.*dq2.*m_p.*t2.*t4.*t13.*t14.*(1.0./2.0);
t37 = m_p.*t7.*t12.*t25.*(1.0./6.0);
t38 = m_p.*t12.*t21.*t25.*(1.0./6.0);
t39 = a.*b.*m_p.*t7.*t15.*(1.0./4.0);
t40 = a.*b.*m_p.*t15.*t21.*(1.0./4.0);
t41 = dq2.*dq3.*m_p.*t4.*t12.*t18.*(2.0./3.0);
t42 = m_p.*t2.*t7.*t11.*t12.*t13.*(1.0./3.0);
t43 = dq1.*dq3.*m_p.*t2.*t12.*t13.*t18.*(2.0./3.0);
t44 = m_p.*t7.*t8.*t9.*t12.*t15.*t20.*(1.0./3.0);
t45 = a.*b.*m_p.*t7.*t8.*t9.*t15.*(1.0./4.0);
t46 = cos(q1);
t47 = sin(q1);
B_ = [-tau1+Iyz.*t3+Iyy.*t3.*t6.*(1.0./2.0)-Iyz.*t3.*t9.*2.0-Iyz.*t7.*t8-Izz.*t3.*t6.*(1.0./2.0)+Iyz.*t7.*t8.*t9.*2.0-m_p.*t3.*t6.*t10.*(1.0./4.0)-m_p.*t3.*t6.*t12.*(1.0./3.0)+m_p.*t3.*t6.*t16.*(1.0./1.2e1)-Ixx.*dq1.*dq2.*t2+Iyy.*dq1.*dq2.*t2-Izz.*dq1.*dq2.*t2+Ixy.*dq1.*dq2.*t11.*t13.*2.0+Ixz.*dq1.*dq2.*t4.*t11.*2.0-Iyy.*dq1.*dq2.*t2.*t9.*2.0+Izz.*dq1.*dq2.*t2.*t9.*2.0-Ixy.*t2.*t4.*t7.*t11+Ixz.*t2.*t7.*t11.*t13-Iyy.*t4.*t7.*t8.*t13+Izz.*t4.*t7.*t8.*t13-dq1.*dq2.*m_p.*t2.*t10-dq1.*dq2.*m_p.*t2.*t12.*(4.0./3.0)+dq2.*dq4.*m_p.*t4.*t12.*(2.0./3.0)-dq2.*dq5.*m_p.*t4.*t12.*(2.0./3.0)+dq3.*dq4.*m_p.*t12.*t23.*(1.0./3.0)+dq3.*dq5.*m_p.*t12.*t25.*(1.0./3.0)+m_p.*t4.*t7.*t8.*t10.*t13.*(1.0./2.0)+m_p.*t4.*t7.*t8.*t12.*t13.*(2.0./3.0)-m_p.*t4.*t7.*t8.*t13.*t16.*(1.0./6.0)+m_p.*t3.*t4.*t12.*t13.*t17.*(1.0./3.0)+m_p.*t3.*t4.*t12.*t13.*t18.*(1.0./3.0)+a.*b.*dq3.*dq4.*m_p.*t14.*(1.0./2.0)+a.*b.*dq3.*dq5.*m_p.*t15.*(1.0./2.0)-Iyz.*dq1.*dq2.*t2.*t4.*t13.*4.0+dq1.*dq2.*m_p.*t2.*t9.*t10+dq1.*dq2.*m_p.*t2.*t9.*t12.*(4.0./3.0)-dq1.*dq2.*m_p.*t2.*t9.*t16.*(1.0./3.0)+dq1.*dq4.*m_p.*t2.*t12.*t13.*(2.0./3.0)-dq1.*dq5.*m_p.*t2.*t12.*t13.*(2.0./3.0)+dq1.*dq2.*m_p.*t2.*t12.*t17.*(2.0./3.0)+dq1.*dq2.*m_p.*t2.*t12.*t18.*(2.0./3.0)-dq2.*dq4.*m_p.*t4.*t12.*t17.*(2.0./3.0)+dq2.*dq5.*m_p.*t4.*t12.*t18.*(2.0./3.0)-a.*b.*dq1.*dq2.*m_p.*t2.*t19-a.*b.*dq1.*dq2.*m_p.*t2.*t20+a.*b.*dq2.*dq4.*m_p.*t4.*t19.*(1.0./2.0)-a.*b.*dq1.*dq4.*m_p.*t11.*t14.*(1.0./2.0)-a.*b.*dq2.*dq5.*m_p.*t4.*t20.*(1.0./2.0)-a.*b.*dq1.*dq5.*m_p.*t11.*t15.*(1.0./2.0)-a.*b.*m_p.*t3.*t4.*t13.*t19.*(1.0./2.0)-a.*b.*m_p.*t3.*t4.*t13.*t20.*(1.0./2.0)-dq1.*dq2.*m_p.*t2.*t9.*t12.*t17.*(2.0./3.0)-dq1.*dq2.*m_p.*t2.*t9.*t12.*t18.*(2.0./3.0)-dq1.*dq4.*m_p.*t2.*t12.*t13.*t17.*(2.0./3.0)+dq1.*dq5.*m_p.*t2.*t12.*t13.*t18.*(2.0./3.0)-dq1.*dq2.*m_p.*t11.*t12.*t13.*t19.*(1.0./2.0)+dq1.*dq2.*m_p.*t11.*t12.*t13.*t20.*(1.0./2.0)-dq1.*dq4.*m_p.*t11.*t12.*t14.*t19.*(2.0./3.0)-dq1.*dq5.*m_p.*t11.*t12.*t15.*t20.*(2.0./3.0)+m_p.*t2.*t4.*t7.*t11.*t12.*t19.*(1.0./4.0)-m_p.*t2.*t4.*t7.*t11.*t12.*t20.*(1.0./4.0)-m_p.*t4.*t7.*t8.*t12.*t13.*t17.*(1.0./3.0)-m_p.*t4.*t7.*t8.*t12.*t13.*t18.*(1.0./3.0)+a.*b.*dq1.*dq2.*m_p.*t2.*t9.*t19+a.*b.*dq1.*dq2.*m_p.*t2.*t9.*t20+a.*b.*dq1.*dq4.*m_p.*t2.*t13.*t19.*(1.0./2.0)+a.*b.*dq1.*dq2.*m_p.*t11.*t13.*t14.*(1.0./2.0)-a.*b.*dq1.*dq5.*m_p.*t2.*t13.*t20.*(1.0./2.0)-a.*b.*dq1.*dq2.*m_p.*t11.*t13.*t15.*(1.0./2.0)-a.*b.*m_p.*t2.*t4.*t7.*t11.*t14.*(1.0./4.0)+a.*b.*m_p.*t2.*t4.*t7.*t11.*t15.*(1.0./4.0)+a.*b.*m_p.*t4.*t7.*t8.*t13.*t19.*(1.0./2.0)+a.*b.*m_p.*t4.*t7.*t8.*t13.*t20.*(1.0./2.0)+dq1.*dq2.*m_p.*t11.*t12.*t13.*t14.*t19.*(2.0./3.0)-dq1.*dq2.*m_p.*t11.*t12.*t13.*t15.*t20.*(2.0./3.0)-m_p.*t2.*t4.*t7.*t11.*t12.*t14.*t19.*(1.0./3.0)+m_p.*t2.*t4.*t7.*t11.*t12.*t15.*t20.*(1.0./3.0);-tau2-Ixz.*t3+Ixz.*t7+Ixz.*t21+Ixy.*t3.*t6.*(1.0./2.0)+Ixz.*t3.*t9-Ixz.*t7.*t8-Ixz.*t7.*t8.*t9-Ixx.*dq2.*dq3.*t13+Ixy.*dq1.*dq2.*t2.*2.0-Ixz.*dq1.*dq3.*t11.*2.0-Iyy.*dq2.*dq3.*t13+Iyz.*dq2.*dq3.*t4.*2.0+Izz.*dq2.*dq3.*t13+Ixx.*dq1.*dq3.*t2.*t4+Ixx.*dq1.*dq2.*t11.*t13-Ixy.*dq1.*dq2.*t2.*t9.*2.0+Iyy.*dq1.*dq3.*t2.*t4-Iyy.*dq1.*dq2.*t11.*t13+Iyz.*dq1.*dq3.*t2.*t13.*2.0-Izz.*dq1.*dq3.*t2.*t4-Izz.*dq1.*dq2.*t11.*t13-Ixx.*t2.*t4.*t7.*t11-Ixy.*t4.*t7.*t8.*t13-Iyz.*t2.*t7.*t11.*t13+Izz.*t2.*t4.*t7.*t11+dq3.*dq4.*m_p.*t12.*t14.*(1.0./2.0)-dq2.*dq3.*m_p.*t13.*t16.*(1.0./3.0)-dq3.*dq5.*m_p.*t12.*t15.*(1.0./2.0)-dq3.*dq4.*m_p.*t12.*t17.*(2.0./3.0)+dq3.*dq5.*m_p.*t12.*t18.*(2.0./3.0)+m_p.*t2.*t4.*t7.*t11.*t12.*(1.0./2.0)-m_p.*t2.*t4.*t7.*t11.*t16.*(1.0./6.0)-m_p.*t3.*t4.*t12.*t13.*t19.*(1.0./4.0)+m_p.*t3.*t4.*t12.*t13.*t20.*(1.0./4.0)+Ixz.*dq1.*dq2.*t2.*t4.*t13.*2.0+dq1.*dq3.*m_p.*t2.*t4.*t16.*(1.0./3.0)-dq1.*dq2.*m_p.*t2.*t12.*t19.*(1.0./2.0)+dq1.*dq2.*m_p.*t2.*t12.*t20.*(1.0./2.0)-dq1.*dq2.*m_p.*t11.*t12.*t13+dq2.*dq4.*m_p.*t4.*t12.*t19.*(1.0./2.0)-dq1.*dq4.*m_p.*t11.*t12.*t14.*(1.0./2.0)+dq2.*dq5.*m_p.*t4.*t12.*t20.*(1.0./2.0)+dq1.*dq5.*m_p.*t11.*t12.*t15.*(1.0./2.0)+dq1.*dq4.*m_p.*t11.*t12.*t17.*(2.0./3.0)-dq1.*dq5.*m_p.*t11.*t12.*t18.*(2.0./3.0)+a.*b.*dq1.*dq2.*m_p.*t2.*t14.*(1.0./2.0)-a.*b.*dq1.*dq2.*m_p.*t2.*t15.*(1.0./2.0)+a.*b.*m_p.*t3.*t4.*t13.*t14.*(1.0./4.0)-a.*b.*m_p.*t3.*t4.*t13.*t15.*(1.0./4.0)+dq1.*dq2.*m_p.*t2.*t9.*t12.*t19.*(1.0./2.0)-dq1.*dq2.*m_p.*t2.*t9.*t12.*t20.*(1.0./2.0)+dq1.*dq2.*m_p.*t2.*t12.*t14.*t19.*(2.0./3.0)+dq1.*dq4.*m_p.*t2.*t12.*t13.*t19.*(1.0./2.0)-dq1.*dq2.*m_p.*t2.*t12.*t15.*t20.*(2.0./3.0)+dq1.*dq2.*m_p.*t11.*t12.*t13.*t14+dq1.*dq5.*m_p.*t2.*t12.*t13.*t20.*(1.0./2.0)+dq1.*dq2.*m_p.*t11.*t12.*t13.*t15-dq2.*dq4.*m_p.*t4.*t12.*t14.*t19.*(2.0./3.0)-dq1.*dq2.*m_p.*t11.*t12.*t13.*t17.*(2.0./3.0)-dq1.*dq2.*m_p.*t11.*t12.*t13.*t18.*(2.0./3.0)-dq2.*dq5.*m_p.*t4.*t12.*t15.*t20.*(2.0./3.0)-m_p.*t2.*t4.*t7.*t11.*t12.*t14.*(1.0./2.0)-m_p.*t2.*t4.*t7.*t11.*t12.*t15.*(1.0./2.0)+m_p.*t2.*t4.*t7.*t11.*t12.*t17.*(1.0./3.0)+m_p.*t2.*t4.*t7.*t11.*t12.*t18.*(1.0./3.0)+m_p.*t4.*t7.*t8.*t12.*t13.*t19.*(1.0./4.0)-m_p.*t4.*t7.*t8.*t12.*t13.*t20.*(1.0./4.0)+m_p.*t3.*t4.*t12.*t13.*t14.*t19.*(1.0./3.0)-m_p.*t3.*t4.*t12.*t13.*t15.*t20.*(1.0./3.0)-a.*b.*dq1.*dq2.*m_p.*t2.*t9.*t14.*(1.0./2.0)+a.*b.*dq1.*dq2.*m_p.*t2.*t9.*t15.*(1.0./2.0)-a.*b.*m_p.*t4.*t7.*t8.*t13.*t14.*(1.0./4.0)+a.*b.*m_p.*t4.*t7.*t8.*t13.*t15.*(1.0./4.0)-dq1.*dq2.*m_p.*t2.*t9.*t12.*t14.*t19.*(2.0./3.0)+dq1.*dq2.*m_p.*t2.*t9.*t12.*t15.*t20.*(2.0./3.0)-dq1.*dq4.*m_p.*t2.*t12.*t13.*t14.*t19.*(2.0./3.0)-dq1.*dq5.*m_p.*t2.*t12.*t13.*t15.*t20.*(2.0./3.0)-m_p.*t4.*t7.*t8.*t12.*t13.*t14.*t19.*(1.0./3.0)+m_p.*t4.*t7.*t8.*t12.*t13.*t15.*t20.*(1.0./3.0);t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t43+t44+t45-tau3-tau_p1-tau_p2-Ixy.*t7-Ixy.*t21+Ixy.*t3.*t9+Ixy.*t7.*t8.*2.0-Ixz.*t3.*t6.*(1.0./2.0)-Ixy.*t7.*t8.*t9+m_p.*t7.*t12.*t19.*(1.0./4.0)-m_p.*t7.*t12.*t20.*(1.0./4.0)-m_p.*t7.*t12.*t23.*(1.0./6.0)+m_p.*t12.*t19.*t21.*(1.0./4.0)-m_p.*t12.*t20.*t21.*(1.0./4.0)-m_p.*t12.*t21.*t23.*(1.0./6.0)+m_p.*t12.*t19.*t26.*(1.0./4.0)-m_p.*t12.*t20.*t27.*(1.0./4.0)-Ixx.*dq2.*dq3.*t4+Ixy.*dq1.*dq3.*t11.*2.0+Iyy.*dq2.*dq3.*t4+Iyz.*dq2.*dq3.*t13.*2.0-Izz.*dq2.*dq3.*t4+Ixx.*dq1.*dq2.*t4.*t11-Ixx.*dq1.*dq3.*t2.*t13+Ixz.*dq1.*dq2.*t2.*t9.*2.0-Iyy.*dq1.*dq2.*t4.*t11+Iyy.*dq1.*dq3.*t2.*t13-Iyz.*dq1.*dq3.*t2.*t4.*2.0-Izz.*dq1.*dq2.*t4.*t11-Izz.*dq1.*dq3.*t2.*t13-a.*b.*m_p.*t7.*t14.*(1.0./4.0)-a.*b.*m_p.*t14.*t21.*(1.0./4.0)+a.*b.*m_p.*t14.*t26.*(1.0./4.0)-a.*b.*m_p.*t15.*t27.*(1.0./4.0)+Ixx.*t2.*t7.*t11.*t13+Ixz.*t4.*t7.*t8.*t13-Iyy.*t2.*t7.*t11.*t13+Iyz.*t2.*t4.*t7.*t11-dq2.*dq3.*m_p.*t4.*t10-dq2.*dq3.*m_p.*t4.*t12.*(4.0./3.0)-m_p.*t3.*t9.*t12.*t19.*(1.0./4.0)+m_p.*t3.*t9.*t12.*t20.*(1.0./4.0)-m_p.*t7.*t8.*t12.*t19.*(1.0./2.0)+m_p.*t7.*t8.*t12.*t20.*(1.0./2.0)+m_p.*t2.*t7.*t10.*t11.*t13.*(1.0./2.0)+m_p.*t2.*t7.*t11.*t12.*t13.*(1.0./6.0)+m_p.*t7.*t8.*t9.*t12.*t19.*(1.0./4.0)-m_p.*t7.*t8.*t9.*t12.*t20.*(1.0./4.0)-m_p.*t3.*t9.*t12.*t15.*t20.*(1.0./3.0)-m_p.*t7.*t8.*t12.*t15.*t20.*(2.0./3.0)+Ixy.*dq1.*dq2.*t2.*t4.*t13.*2.0-a.*b.*m_p.*t3.*t9.*t15.*(1.0./4.0)+a.*b.*m_p.*t7.*t8.*t14.*(1.0./2.0)-a.*b.*m_p.*t7.*t8.*t15.*(1.0./2.0)-dq1.*dq3.*m_p.*t2.*t10.*t13-dq1.*dq2.*m_p.*t4.*t11.*t12-dq1.*dq3.*m_p.*t2.*t12.*t13.*(4.0./3.0)-dq1.*dq3.*m_p.*t11.*t12.*t19.*(1.0./2.0)+dq1.*dq3.*m_p.*t11.*t12.*t20.*(1.0./2.0)-dq2.*dq4.*m_p.*t12.*t13.*t19.*(1.0./2.0)-dq2.*dq5.*m_p.*t12.*t13.*t20.*(1.0./2.0)-a.*b.*dq2.*dq3.*m_p.*t4.*t19-a.*b.*dq2.*dq3.*m_p.*t4.*t20-a.*b.*dq1.*dq3.*m_p.*t11.*t15.*(1.0./2.0)-a.*b.*dq2.*dq4.*m_p.*t13.*t14.*(1.0./2.0)-a.*b.*dq2.*dq5.*m_p.*t13.*t15.*(1.0./2.0)-a.*b.*m_p.*t7.*t8.*t9.*t14.*(1.0./4.0)+dq1.*dq4.*m_p.*t2.*t4.*t12.*t19.*(1.0./2.0)+dq1.*dq2.*m_p.*t4.*t11.*t12.*t14+dq1.*dq5.*m_p.*t2.*t4.*t12.*t20.*(1.0./2.0)+dq1.*dq2.*m_p.*t4.*t11.*t12.*t15-dq1.*dq2.*m_p.*t4.*t11.*t12.*t17.*(2.0./3.0)-dq1.*dq2.*m_p.*t4.*t11.*t12.*t18.*(2.0./3.0)-dq1.*dq3.*m_p.*t11.*t12.*t15.*t20.*(2.0./3.0)+m_p.*t2.*t7.*t11.*t12.*t13.*t14.*(1.0./2.0)+m_p.*t2.*t7.*t11.*t12.*t13.*t15.*(1.0./2.0)-m_p.*t2.*t7.*t11.*t12.*t13.*t17.*(2.0./3.0)-m_p.*t2.*t7.*t11.*t12.*t13.*t18.*(2.0./3.0)-m_p.*t7.*t8.*t9.*t12.*t14.*t19.*(1.0./3.0)+a.*b.*dq1.*dq4.*m_p.*t2.*t4.*t14.*(1.0./2.0)+a.*b.*dq1.*dq5.*m_p.*t2.*t4.*t15.*(1.0./2.0)-a.*b.*dq1.*dq3.*m_p.*t2.*t13.*t19-a.*b.*dq1.*dq3.*m_p.*t2.*t13.*t20+a.*b.*m_p.*t2.*t7.*t11.*t13.*t19.*(1.0./2.0)+a.*b.*m_p.*t2.*t7.*t11.*t13.*t20.*(1.0./2.0)-dq1.*dq2.*m_p.*t2.*t4.*t12.*t13.*t19.*(1.0./2.0)+dq1.*dq2.*m_p.*t2.*t4.*t12.*t13.*t20.*(1.0./2.0)-a.*b.*dq1.*dq2.*m_p.*t2.*t4.*t13.*t15.*(1.0./2.0)-dq1.*dq2.*m_p.*t2.*t4.*t12.*t13.*t15.*t20.*(2.0./3.0);t28+t29+t30+t31+t32+t33+t34+t35+t36+t42-tau_p1-m_p.*t3.*t12.*t19.*(1.0./4.0)-m_p.*t7.*t12.*t23.*(1.0./6.0)-m_p.*t12.*t21.*t23.*(1.0./6.0)-a.*b.*m_p.*t3.*t14.*(1.0./4.0)-a.*b.*m_p.*t7.*t14.*(1.0./4.0)-a.*b.*m_p.*t14.*t21.*(1.0./4.0)-dq2.*dq3.*m_p.*t4.*t12.*(2.0./3.0)-m_p.*t7.*t8.*t12.*t19.*(1.0./4.0)+a.*b.*m_p.*t7.*t8.*t14.*(1.0./4.0)-dq1.*dq3.*m_p.*t2.*t12.*t13.*(2.0./3.0)-a.*b.*dq2.*dq3.*m_p.*t4.*t19.*(1.0./2.0)-a.*b.*m_p.*t7.*t8.*t9.*t14.*(1.0./4.0)+dq1.*dq2.*m_p.*t4.*t11.*t12.*t14.*(1.0./2.0)-dq1.*dq2.*m_p.*t4.*t11.*t12.*t17.*(2.0./3.0)+m_p.*t2.*t7.*t11.*t12.*t13.*t14.*(1.0./4.0)-m_p.*t2.*t7.*t11.*t12.*t13.*t17.*(2.0./3.0)-m_p.*t7.*t8.*t9.*t12.*t14.*t19.*(1.0./3.0)-a.*b.*dq1.*dq3.*m_p.*t2.*t13.*t19.*(1.0./2.0)+a.*b.*m_p.*t2.*t7.*t11.*t13.*t19.*(1.0./4.0);-t37-t38-t39-t40-t41-t42-t43-t44-t45+tau_p2-m_p.*t3.*t12.*t20.*(1.0./4.0)-a.*b.*m_p.*t3.*t15.*(1.0./4.0)+dq2.*dq3.*m_p.*t4.*t12.*(2.0./3.0)-m_p.*t7.*t8.*t12.*t20.*(1.0./4.0)+m_p.*t3.*t9.*t12.*t15.*t20.*(1.0./3.0)+m_p.*t7.*t8.*t12.*t15.*t20.*(2.0./3.0)+a.*b.*m_p.*t3.*t9.*t15.*(1.0./4.0)+a.*b.*m_p.*t7.*t8.*t15.*(1.0./4.0)+dq1.*dq3.*m_p.*t2.*t12.*t13.*(2.0./3.0)+a.*b.*dq2.*dq3.*m_p.*t4.*t20.*(1.0./2.0)+a.*b.*dq1.*dq3.*m_p.*t11.*t15.*(1.0./2.0)-dq1.*dq2.*m_p.*t4.*t11.*t12.*t15.*(1.0./2.0)+dq1.*dq2.*m_p.*t4.*t11.*t12.*t18.*(2.0./3.0)+dq1.*dq3.*m_p.*t11.*t12.*t15.*t20.*(2.0./3.0)-m_p.*t2.*t7.*t11.*t12.*t13.*t15.*(1.0./4.0)+m_p.*t2.*t7.*t11.*t12.*t13.*t18.*(2.0./3.0)+a.*b.*dq1.*dq3.*m_p.*t2.*t13.*t20.*(1.0./2.0)-a.*b.*m_p.*t2.*t7.*t11.*t13.*t20.*(1.0./4.0)+a.*b.*dq1.*dq2.*m_p.*t2.*t4.*t13.*t15.*(1.0./2.0)+dq1.*dq2.*m_p.*t2.*t4.*t12.*t13.*t15.*t20.*(2.0./3.0);b.*m_p.*(t2.*t3.*t46.*2.0+t2.*t7.*t46.*2.0-dq1.*dq2.*t11.*t47.*4.0-t2.*t3.*t14.*t46-t2.*t3.*t15.*t46-t2.*t7.*t14.*t46-t2.*t7.*t15.*t46+t4.*t7.*t19.*t47-t4.*t7.*t20.*t47-t2.*t14.*t26.*t46-t2.*t15.*t27.*t46+t4.*t19.*t21.*t47-t4.*t20.*t21.*t47+t4.*t19.*t26.*t47-t4.*t20.*t27.*t47-dq1.*dq4.*t4.*t14.*t46.*2.0+dq1.*dq5.*t4.*t15.*t46.*2.0+dq1.*dq4.*t2.*t19.*t47.*2.0+dq1.*dq2.*t11.*t14.*t47.*2.0+dq1.*dq5.*t2.*t20.*t47.*2.0+dq1.*dq2.*t11.*t15.*t47.*2.0+dq3.*dq4.*t13.*t14.*t47.*2.0+dq1.*dq3.*t13.*t19.*t46.*2.0+dq2.*dq4.*t11.*t19.*t46.*2.0-dq1.*dq3.*t13.*t20.*t46.*2.0-dq3.*dq5.*t13.*t15.*t47.*2.0+dq2.*dq5.*t11.*t20.*t46.*2.0-t3.*t11.*t13.*t19.*t46+t3.*t11.*t13.*t20.*t46-t7.*t11.*t13.*t19.*t46+t7.*t11.*t13.*t20.*t46-t11.*t13.*t19.*t21.*t46+t11.*t13.*t20.*t21.*t46-t11.*t13.*t19.*t26.*t46+t11.*t13.*t20.*t27.*t46+dq2.*dq3.*t2.*t4.*t19.*t46.*2.0-dq2.*dq3.*t2.*t4.*t20.*t46.*2.0+dq2.*dq4.*t2.*t13.*t14.*t46.*2.0+dq3.*dq4.*t4.*t11.*t14.*t46.*2.0-dq2.*dq5.*t2.*t13.*t15.*t46.*2.0-dq1.*dq2.*t2.*t13.*t19.*t47.*2.0-dq3.*dq5.*t4.*t11.*t15.*t46.*2.0+dq1.*dq2.*t2.*t13.*t20.*t47.*2.0-dq1.*dq3.*t4.*t11.*t19.*t47.*2.0+dq1.*dq3.*t4.*t11.*t20.*t47.*2.0-dq1.*dq4.*t11.*t13.*t14.*t47.*2.0+dq1.*dq5.*t11.*t13.*t15.*t47.*2.0).*(-1.0./2.0);b.*m_p.*(t2.*t3.*t47.*-2.0-t2.*t7.*t47.*2.0-dq1.*dq2.*t11.*t46.*4.0+t2.*t3.*t14.*t47+t2.*t3.*t15.*t47+t2.*t7.*t14.*t47+t2.*t7.*t15.*t47+t4.*t7.*t19.*t46-t4.*t7.*t20.*t46+t2.*t14.*t26.*t47+t4.*t19.*t21.*t46+t2.*t15.*t27.*t47-t4.*t20.*t21.*t46+t4.*t19.*t26.*t46-t4.*t20.*t27.*t46+dq1.*dq4.*t4.*t14.*t47.*2.0+dq1.*dq4.*t2.*t19.*t46.*2.0-dq1.*dq5.*t4.*t15.*t47.*2.0+dq1.*dq2.*t11.*t14.*t46.*2.0+dq1.*dq5.*t2.*t20.*t46.*2.0+dq1.*dq2.*t11.*t15.*t46.*2.0+dq3.*dq4.*t13.*t14.*t46.*2.0-dq3.*dq5.*t13.*t15.*t46.*2.0-dq1.*dq3.*t13.*t19.*t47.*2.0-dq2.*dq4.*t11.*t19.*t47.*2.0+dq1.*dq3.*t13.*t20.*t47.*2.0-dq2.*dq5.*t11.*t20.*t47.*2.0+t3.*t11.*t13.*t19.*t47-t3.*t11.*t13.*t20.*t47+t7.*t11.*t13.*t19.*t47-t7.*t11.*t13.*t20.*t47+t11.*t13.*t19.*t21.*t47-t11.*t13.*t20.*t21.*t47+t11.*t13.*t19.*t26.*t47-t11.*t13.*t20.*t27.*t47-dq2.*dq3.*t2.*t4.*t19.*t47.*2.0+dq2.*dq3.*t2.*t4.*t20.*t47.*2.0-dq2.*dq4.*t2.*t13.*t14.*t47.*2.0-dq1.*dq2.*t2.*t13.*t19.*t46.*2.0-dq3.*dq4.*t4.*t11.*t14.*t47.*2.0+dq1.*dq2.*t2.*t13.*t20.*t46.*2.0-dq1.*dq3.*t4.*t11.*t19.*t46.*2.0+dq2.*dq5.*t2.*t13.*t15.*t47.*2.0+dq1.*dq3.*t4.*t11.*t20.*t46.*2.0+dq3.*dq5.*t4.*t11.*t15.*t47.*2.0-dq1.*dq4.*t11.*t13.*t14.*t46.*2.0+dq1.*dq5.*t11.*t13.*t15.*t46.*2.0).*(1.0./2.0);b.*m_p.*(t3.*t11.*-2.0+t3.*t11.*t14+t3.*t11.*t15+t11.*t14.*t26+t11.*t15.*t27+dq2.*dq4.*t2.*t19.*2.0+dq2.*dq5.*t2.*t20.*2.0-t2.*t3.*t13.*t19+t2.*t3.*t13.*t20-t2.*t13.*t19.*t21+t2.*t13.*t20.*t21-t2.*t13.*t19.*t26+t2.*t13.*t20.*t27+dq3.*dq4.*t2.*t4.*t14.*2.0-dq3.*dq5.*t2.*t4.*t15.*2.0-dq2.*dq3.*t4.*t11.*t19.*2.0+dq2.*dq3.*t4.*t11.*t20.*2.0-dq2.*dq4.*t11.*t13.*t14.*2.0+dq2.*dq5.*t11.*t13.*t15.*2.0).*(-1.0./2.0)];