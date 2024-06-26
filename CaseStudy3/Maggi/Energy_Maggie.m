function T = Energy_Maggie(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,dq1,dq2,dq3,dq4,dq5,dq6,dq7,l,m,m_sat,q1,q2,q3,q4)
%ENERGY_MAGGIE
%    T = ENERGY_MAGGIE(IXX,IXY,IXZ,IYY,IYZ,IZZ,A,DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,L,M,M_SAT,Q1,Q2,Q3,Q4)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    07-Oct-2020 06:47:13

t2 = cos(q2);
t4 = sin(q1);
t5 = cos(q1);
t6 = sin(q2);
t39 = dq7.*t6;
t42 = dq5.*t2.*t5;
t43 = dq6.*t2.*t4;
t3 = -t39+t42+t43;
t7 = cos(q3);
t8 = sin(q3);
t11 = t4.*t7;
t12 = t5.*t6.*t8;
t13 = t11-t12;
t14 = dq5.*t13;
t15 = t5.*t7;
t16 = t4.*t6.*t8;
t17 = t15+t16;
t18 = dq6.*t17;
t19 = dq7.*t2.*t8;
t21 = dq1.*t6;
t22 = dq3-t21;
t40 = l.*(1.0./2.0);
t41 = a+q4+t40;
t9 = t14-t18-t19+t22.*t41;
t28 = t4.*t8;
t29 = t5.*t6.*t7;
t30 = t28+t29;
t31 = dq5.*t30;
t32 = t5.*t8;
t33 = t4.*t6.*t7;
t34 = t32-t33;
t35 = dq6.*t34;
t36 = dq7.*t2.*t7;
t10 = t31-t35+t36;
t20 = -t14+t18+t19;
t23 = dq2.*t7;
t24 = dq1.*t2.*t8;
t25 = t23+t24;
t26 = dq2.*t8;
t38 = dq1.*t2.*t7;
t27 = t26-t38;
t37 = dq4+t31-t35+t36;
t44 = -t39+t42+t43+t25.*t41;
T = (dq2.*t7.*(1.0./2.0)+dq1.*t2.*t8.*(1.0./2.0)).*(-Ixy.*t22+Iyy.*t25+Iyz.*t27)+(dq2.*t8.*(1.0./2.0)-dq1.*t2.*t7.*(1.0./2.0)).*(Ixz.*t22+Iyz.*t25+Izz.*t27)+(dq3.*(1.0./2.0)-dq1.*t6.*(1.0./2.0)).*(Ixx.*t22-Ixy.*t25+Ixz.*t27)+m.*t9.^2.*(1.0./2.0)+m.*t37.^2.*(1.0./2.0)+m.*t44.^2.*(1.0./2.0)+m_sat.*t3.^2.*(1.0./2.0)+m_sat.*t10.^2.*(1.0./2.0)+m_sat.*t20.^2.*(1.0./2.0);
