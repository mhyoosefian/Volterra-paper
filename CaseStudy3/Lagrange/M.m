function M = M(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,l,m,m_sat,q1,q2,q3,q4)
%M
%    M = M(IXX,IXY,IXZ,IYY,IYZ,IZZ,A,L,M,M_SAT,Q1,Q2,Q3,Q4)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    27-Sep-2020 22:01:57

t2 = sin(q2);
t3 = cos(q2);
t7 = l.*(1.0./2.0);
t4 = a+q4+t7;
t5 = cos(q3);
t6 = sin(q3);
t8 = t4.^2;
t9 = Ixy.*t2;
t10 = Iyy.*t3.*t6;
t24 = Iyz.*t3.*t5;
t11 = t9+t10-t24;
t12 = Ixz.*t2;
t13 = Izz.*t3.*t5;
t26 = Iyz.*t3.*t6;
t14 = t12+t13-t26;
t15 = Ixx.*t2;
t16 = Ixz.*t3.*t5;
t17 = Ixy.*t3.*t6;
t18 = cos(q1);
t19 = t3.^2;
t20 = sin(q1);
t21 = Ixy.*t5;
t36 = Ixz.*t6;
t22 = t21-t36;
t23 = t2.*t22.*(1.0./2.0);
t25 = t5.*t11.*(1.0./2.0);
t27 = Iyz.*t5;
t28 = Izz.*t6;
t29 = t27+t28;
t30 = Iyy.*t5;
t31 = Iyz.*t6;
t32 = t30+t31;
t33 = t3.*t6.*t32.*(1.0./2.0);
t34 = m.*t3.*t5.*t6.*t8;
t35 = t23+t25+t33+t34-t6.*t14.*(1.0./2.0)-t3.*t5.*t29.*(1.0./2.0);
t37 = -t15-t16-t17-m.*t2.*t8;
t38 = -t21+t36;
t39 = t5.*t20;
t44 = t2.*t6.*t18;
t40 = t39-t44;
t41 = t5.*t18;
t42 = t2.*t6.*t20;
t43 = t41+t42;
t45 = m.*t4.*t6.*t18.*t19;
t46 = t45-m.*t2.*t4.*t40;
t47 = m.*t3.*t4.*t5.*t18;
t48 = m.*t4.*t40;
t49 = t6.*t20;
t50 = t2.*t5.*t18;
t51 = t49+t50;
t52 = t51.^2;
t53 = t40.^2;
t54 = t18.^2;
t55 = t6.*t18;
t61 = t2.*t5.*t20;
t56 = t55-t61;
t57 = m.*t2.*t4.*t43;
t58 = m.*t4.*t6.*t19.*t20;
t59 = t57+t58;
t60 = m.*t3.*t4.*t5.*t20;
t62 = m.*t18.*t19.*t20;
t63 = m_sat.*t18.*t19.*t20;
t64 = t62+t63-m.*t40.*t43-m.*t51.*t56-m_sat.*t40.*t43-m_sat.*t51.*t56;
t65 = t43.^2;
t66 = t56.^2;
t67 = t20.^2;
t68 = m.*t3.*t5;
t69 = m.*t3.*t5.*t51;
t70 = m_sat.*t3.*t5.*t51;
t71 = t69+t70-m.*t2.*t3.*t18-m.*t3.*t6.*t40-m_sat.*t2.*t3.*t18-m_sat.*t3.*t6.*t40;
t72 = m.*t3.*t6.*t43;
t73 = m_sat.*t3.*t6.*t43;
t74 = t72+t73-m.*t2.*t3.*t20-m.*t3.*t5.*t56-m_sat.*t2.*t3.*t20-m_sat.*t3.*t5.*t56;
t75 = t2.^2;
t76 = t5.^2;
t77 = t6.^2;
M = reshape([t2.*(t15+t16+t17)+m.*t8.*t75+t3.*t6.*t11+t3.*t5.*t14+m.*t8.*t19.*t77,t35,t37,0.0,t46,t59,0.0,t35,t6.*t29+t5.*t32+m.*t8.*t76,t38,0.0,t47,t60,-m.*t2.*t4.*t5,t37,t38,Ixx+m.*t8,0.0,t48,-m.*t4.*t43,-m.*t3.*t4.*t6,0.0,0.0,0.0,m,m.*t51,-m.*t56,t68,t46,t47,t48,m.*(t6.*t20.*2.0+t2.*t5.*t18.*2.0).*(1.0./2.0),m.*t52+m.*t53+m_sat.*t52+m_sat.*t53+m.*t19.*t54+m_sat.*t19.*t54,t64,t71,t59,t60,-m.*t4.*t43,m.*(t6.*t18.*2.0-t2.*t5.*t20.*2.0).*(-1.0./2.0),t64,m.*t65+m.*t66+m_sat.*t65+m_sat.*t66+m.*t19.*t67+m_sat.*t19.*t67,t74,0.0,-m.*t2.*t4.*t5,-m.*t3.*t4.*t6,t68,t71,t74,m.*t75+m_sat.*t75+m.*t19.*t76+m.*t19.*t77+m_sat.*t19.*t76+m_sat.*t19.*t77],[7,7]);
