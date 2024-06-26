function B_ = Bias(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,U5,a,b,c,m_p,m_s,q4,q5,tau1,tau2,tau3,tau_p1,tau_p2)
%BIAS
%    B_ = BIAS(IXX,IXY,IXZ,IYY,IYZ,IZZ,U1,U2,U3,U4,U5,A,B,C,M_P,M_S,Q4,Q5,TAU1,TAU2,TAU3,TAU_P1,TAU_P2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    26-Aug-2020 09:59:22

t2 = U2.^2;
t3 = U3.^2;
t4 = m_p.^2;
t5 = b.^2;
t6 = a.^2;
t7 = c.^2;
t8 = cos(q4);
t9 = cos(q5);
t10 = t8.^2;
t11 = t9.^2;
t12 = q4.*2.0;
t13 = sin(t12);
t14 = q5.*2.0;
t15 = sin(t14);
t16 = sin(q5);
t17 = sin(q4);
t18 = U1.^2;
t19 = m_p.*2.4e1;
t20 = m_s.*1.2e1;
t21 = t19+t20;
t22 = 1.0./t21;
t23 = q4+q5;
t24 = sin(t23);
t25 = q4-q5;
t26 = sin(t25);
t27 = U4.^2;
t28 = U5.^2;
t29 = cos(t12);
t30 = cos(t14);
t31 = m_p.*4.8e1;
t32 = m_s.*2.4e1;
t33 = t31+t32;
t34 = 1.0./t33;
t35 = t2.*t4.*t5.*t26.*3.0;
t36 = cos(t25);
t37 = U1.*U2.*t4.*t5.*t36.*6.0;
t38 = U1.*U2.*m_p.*m_s.*t5.*t8.*6.0;
t39 = U1.*U2.*a.*b.*t4.*t17.*1.2e1;
t40 = U1.*U2.*a.*b.*m_p.*m_s.*t17.*6.0;
t41 = t4.*t5.*t18.*t24.*3.0;
t42 = t2.*t4.*t5.*t24.*3.0;
t43 = t3.*t4.*t5.*t24.*6.0;
t44 = U1.*U2.*m_p.*m_s.*t5.*t9.*6.0;
t45 = U1.*U2.*a.*b.*t4.*t16.*1.2e1;
t46 = U1.*U2.*a.*b.*m_p.*m_s.*t16.*6.0;
B_ = [t22.*(m_p.*tau1.*-2.4e1-m_s.*tau1.*1.2e1-Iyz.*m_p.*t2.*2.4e1+Iyz.*m_p.*t3.*2.4e1-Iyz.*m_s.*t2.*1.2e1+Iyz.*m_s.*t3.*1.2e1+Ixy.*U1.*U3.*m_p.*2.4e1-Ixz.*U1.*U2.*m_p.*2.4e1+Ixy.*U1.*U3.*m_s.*1.2e1-Ixz.*U1.*U2.*m_s.*1.2e1-Iyy.*U2.*U3.*m_p.*2.4e1-Iyy.*U2.*U3.*m_s.*1.2e1+Izz.*U2.*U3.*m_p.*2.4e1+Izz.*U2.*U3.*m_s.*1.2e1+U2.*U3.*t4.*t5.*1.0e1+U2.*U3.*t4.*t6.*1.2e1+U2.*U4.*t4.*t5.*1.0e1-U2.*U3.*t4.*t7.*4.0-U2.*U5.*t4.*t5.*1.0e1+U2.*U3.*m_p.*m_s.*t5.*8.0+U2.*U3.*m_p.*m_s.*t6.*6.0+U2.*U4.*m_p.*m_s.*t5.*8.0-U2.*U3.*m_p.*m_s.*t7.*2.0-U2.*U5.*m_p.*m_s.*t5.*8.0-U2.*U3.*t4.*t5.*t10.*5.0-U2.*U3.*t4.*t5.*t11.*5.0-U2.*U4.*t4.*t5.*t10.*1.0e1+U1.*U3.*t4.*t5.*t13.*(5.0./2.0)+U1.*U4.*t4.*t5.*t13.*5.0+U2.*U5.*t4.*t5.*t11.*1.0e1-U1.*U3.*t4.*t5.*t15.*(5.0./2.0)+U1.*U5.*t4.*t5.*t15.*5.0+U1.*U3.*a.*b.*t4.*t8.*6.0-U1.*U3.*a.*b.*t4.*t9.*6.0+U1.*U4.*a.*b.*t4.*t8.*1.2e1+U1.*U5.*a.*b.*t4.*t9.*1.2e1+U2.*U3.*a.*b.*t4.*t16.*1.2e1+U2.*U3.*a.*b.*t4.*t17.*1.2e1+U2.*U4.*a.*b.*t4.*t17.*1.2e1-U2.*U5.*a.*b.*t4.*t16.*1.2e1-U2.*U3.*m_p.*m_s.*t5.*t10.*4.0-U2.*U3.*m_p.*m_s.*t5.*t11.*4.0-U2.*U4.*m_p.*m_s.*t5.*t10.*8.0+U1.*U3.*m_p.*m_s.*t5.*t13.*2.0+U1.*U4.*m_p.*m_s.*t5.*t13.*4.0+U2.*U5.*m_p.*m_s.*t5.*t11.*8.0-U1.*U3.*m_p.*m_s.*t5.*t15.*2.0+U1.*U3.*m_p.*m_s.*t5.*t16.*3.0-U1.*U3.*m_p.*m_s.*t5.*t17.*3.0+U1.*U5.*m_p.*m_s.*t5.*t15.*4.0+U1.*U3.*t4.*t5.*t8.*t16.*3.0+U1.*U4.*t4.*t5.*t8.*t16.*6.0-U1.*U3.*t4.*t5.*t9.*t17.*3.0+U1.*U5.*t4.*t5.*t9.*t17.*6.0+U2.*U3.*t4.*t5.*t16.*t17.*6.0+U2.*U4.*t4.*t5.*t16.*t17.*6.0-U2.*U5.*t4.*t5.*t16.*t17.*6.0+U1.*U3.*a.*b.*m_p.*m_s.*t8.*3.0-U1.*U3.*a.*b.*m_p.*m_s.*t9.*3.0+U1.*U4.*a.*b.*m_p.*m_s.*t8.*6.0+U1.*U5.*a.*b.*m_p.*m_s.*t9.*6.0+U2.*U3.*a.*b.*m_p.*m_s.*t16.*6.0+U2.*U3.*a.*b.*m_p.*m_s.*t17.*6.0+U2.*U4.*a.*b.*m_p.*m_s.*t17.*6.0-U2.*U5.*a.*b.*m_p.*m_s.*t16.*6.0);-t34.*(m_p.*tau2.*4.8e1+m_s.*tau2.*2.4e1+Ixz.*m_p.*t3.*4.8e1-Ixz.*m_p.*t18.*4.8e1+Ixz.*m_s.*t3.*2.4e1-Ixz.*m_s.*t18.*2.4e1-Ixx.*U1.*U3.*m_p.*4.8e1+Ixy.*U2.*U3.*m_p.*4.8e1-Ixx.*U1.*U3.*m_s.*2.4e1+Ixy.*U2.*U3.*m_s.*2.4e1-Iyz.*U1.*U2.*m_p.*4.8e1-Iyz.*U1.*U2.*m_s.*2.4e1+Izz.*U1.*U3.*m_p.*4.8e1+Izz.*U1.*U3.*m_s.*2.4e1-U1.*U3.*t4.*t7.*8.0+U1.*U3.*m_p.*m_s.*t5.*1.2e1-U1.*U3.*m_p.*m_s.*t7.*4.0+U1.*U3.*t4.*t5.*t10.*1.0e1+U1.*U3.*t4.*t5.*t11.*1.0e1+U1.*U4.*t4.*t5.*t10.*2.0e1-U1.*U5.*t4.*t5.*t11.*2.0e1+U2.*U3.*t4.*t5.*t13.*5.0+U2.*U4.*t4.*t5.*t13.*1.0e1-U2.*U3.*t4.*t5.*t15.*5.0+U2.*U5.*t4.*t5.*t15.*1.0e1+U2.*U3.*a.*b.*t4.*t8.*1.2e1-U2.*U3.*a.*b.*t4.*t9.*1.2e1-U1.*U3.*m_p.*m_s.*t5.*t8.*1.2e1-U1.*U3.*m_p.*m_s.*t5.*t9.*1.2e1-U1.*U4.*m_p.*m_s.*t5.*t8.*1.2e1+U1.*U3.*m_p.*m_s.*t5.*t10.*8.0+U1.*U3.*m_p.*m_s.*t5.*t11.*8.0+U1.*U4.*m_p.*m_s.*t5.*t10.*1.6e1+U1.*U5.*m_p.*m_s.*t5.*t9.*1.2e1-U1.*U5.*m_p.*m_s.*t5.*t11.*1.6e1+U2.*U3.*m_p.*m_s.*t5.*t13.*4.0+U2.*U4.*m_p.*m_s.*t5.*t13.*8.0-U2.*U3.*m_p.*m_s.*t5.*t15.*4.0+U2.*U3.*m_p.*m_s.*t5.*t16.*6.0-U2.*U3.*m_p.*m_s.*t5.*t17.*6.0+U2.*U5.*m_p.*m_s.*t5.*t15.*8.0-U2.*U4.*m_p.*m_s.*t5.*t17.*1.2e1-U2.*U5.*m_p.*m_s.*t5.*t16.*1.2e1-U1.*U3.*t4.*t5.*t8.*t9.*1.2e1-U1.*U4.*t4.*t5.*t8.*t9.*1.2e1+U1.*U5.*t4.*t5.*t8.*t9.*1.2e1+U2.*U3.*t4.*t5.*t8.*t16.*6.0-U2.*U3.*t4.*t5.*t9.*t17.*6.0-U2.*U5.*t4.*t5.*t8.*t16.*1.2e1-U2.*U4.*t4.*t5.*t9.*t17.*1.2e1+U2.*U3.*a.*b.*m_p.*m_s.*t8.*6.0-U2.*U3.*a.*b.*m_p.*m_s.*t9.*6.0);-t22.*(t35+t37+t38+t39+t40+t44+t45+t46+m_p.*tau3.*2.4e1+m_p.*tau_p1.*2.4e1+m_p.*tau_p2.*2.4e1+m_s.*tau3.*1.2e1+m_s.*tau_p1.*1.2e1+m_s.*tau_p2.*1.2e1-Ixy.*m_p.*t2.*2.4e1+Ixy.*m_p.*t18.*2.4e1-Ixy.*m_s.*t2.*1.2e1+Ixy.*m_s.*t18.*1.2e1-t2.*t4.*t5.*t13.*(5.0./2.0)+t2.*t4.*t5.*t15.*(5.0./2.0)+t4.*t5.*t13.*t18.*(5.0./2.0)-t4.*t5.*t15.*t18.*(5.0./2.0)-t4.*t5.*t18.*t26.*3.0-t4.*t5.*t24.*t27.*3.0+t4.*t5.*t24.*t28.*3.0+Ixx.*U1.*U2.*m_p.*2.4e1-Ixz.*U2.*U3.*m_p.*2.4e1+Ixx.*U1.*U2.*m_s.*1.2e1-Ixz.*U2.*U3.*m_s.*1.2e1-Iyy.*U1.*U2.*m_p.*2.4e1+Iyz.*U1.*U3.*m_p.*2.4e1-Iyy.*U1.*U2.*m_s.*1.2e1+Iyz.*U1.*U3.*m_s.*1.2e1+U1.*U2.*t4.*t6.*1.2e1-U1.*U2.*m_p.*m_s.*t5.*6.0+U1.*U2.*m_p.*m_s.*t6.*6.0-U3.*U4.*t4.*t5.*t24.*6.0-U1.*U2.*t4.*t5.*t29.*5.0-U3.*U5.*t4.*t5.*t24.*6.0-U1.*U2.*t4.*t5.*t30.*5.0-a.*b.*t2.*t4.*t8.*6.0+a.*b.*t2.*t4.*t9.*6.0+a.*b.*t4.*t8.*t18.*6.0-a.*b.*t4.*t9.*t18.*6.0-a.*b.*t4.*t8.*t27.*6.0+a.*b.*t4.*t9.*t28.*6.0-m_p.*m_s.*t2.*t5.*t13.*2.0+m_p.*m_s.*t2.*t5.*t15.*2.0-m_p.*m_s.*t2.*t5.*t16.*3.0+m_p.*m_s.*t2.*t5.*t17.*3.0+m_p.*m_s.*t5.*t13.*t18.*2.0-m_p.*m_s.*t5.*t15.*t18.*2.0+m_p.*m_s.*t5.*t16.*t18.*3.0-m_p.*m_s.*t5.*t17.*t18.*3.0+m_p.*m_s.*t5.*t16.*t28.*3.0-m_p.*m_s.*t5.*t17.*t27.*3.0-U3.*U4.*a.*b.*t4.*t8.*1.2e1-U3.*U5.*a.*b.*t4.*t9.*1.2e1-U3.*U4.*m_p.*m_s.*t5.*t17.*6.0-U3.*U5.*m_p.*m_s.*t5.*t16.*6.0-U1.*U2.*m_p.*m_s.*t5.*t29.*4.0-U1.*U2.*m_p.*m_s.*t5.*t30.*4.0-a.*b.*m_p.*m_s.*t2.*t8.*3.0+a.*b.*m_p.*m_s.*t2.*t9.*3.0+a.*b.*m_p.*m_s.*t8.*t18.*3.0-a.*b.*m_p.*m_s.*t9.*t18.*3.0-a.*b.*m_p.*m_s.*t8.*t27.*3.0+a.*b.*m_p.*m_s.*t9.*t28.*3.0-U3.*U4.*a.*b.*m_p.*m_s.*t8.*6.0-U3.*U5.*a.*b.*m_p.*m_s.*t9.*6.0);-t34.*(t35+t37+t38+t39+t40+t41+t42+t43+m_p.*tau_p1.*4.8e1+m_s.*tau_p1.*2.4e1-t2.*t4.*t5.*t13.*5.0+t4.*t5.*t13.*t18.*5.0-t4.*t5.*t18.*t26.*3.0+t4.*t5.*t24.*t28.*6.0-U1.*U2.*t4.*t5.*t29.*1.0e1-U3.*U5.*t4.*t5.*t24.*1.2e1+a.*b.*t3.*t4.*t8.*1.2e1+a.*b.*t4.*t8.*t18.*1.2e1-m_p.*m_s.*t2.*t5.*t13.*4.0+m_p.*m_s.*t2.*t5.*t17.*6.0+m_p.*m_s.*t3.*t5.*t17.*6.0+m_p.*m_s.*t5.*t13.*t18.*4.0-U1.*U2.*m_p.*m_s.*t5.*t29.*8.0+a.*b.*m_p.*m_s.*t3.*t8.*6.0+a.*b.*m_p.*m_s.*t8.*t18.*6.0);-t34.*(-t35-t37+t41+t42+t43-t44-t45-t46-m_p.*tau_p2.*4.8e1-m_s.*tau_p2.*2.4e1-t2.*t4.*t5.*t15.*5.0+t4.*t5.*t15.*t18.*5.0+t4.*t5.*t18.*t26.*3.0+t4.*t5.*t24.*t27.*6.0+U3.*U4.*t4.*t5.*t24.*1.2e1+U1.*U2.*t4.*t5.*t30.*1.0e1+a.*b.*t3.*t4.*t9.*1.2e1+a.*b.*t4.*t9.*t18.*1.2e1-m_p.*m_s.*t2.*t5.*t15.*4.0+m_p.*m_s.*t2.*t5.*t16.*6.0+m_p.*m_s.*t3.*t5.*t16.*6.0+m_p.*m_s.*t5.*t15.*t18.*4.0+U1.*U2.*m_p.*m_s.*t5.*t30.*8.0+a.*b.*m_p.*m_s.*t3.*t9.*6.0+a.*b.*m_p.*m_s.*t9.*t18.*6.0)];
