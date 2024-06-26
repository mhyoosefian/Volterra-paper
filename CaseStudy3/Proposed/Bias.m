function B_ = Bias(F,Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,a,l,m,m_sat,q4,tau1,tau2,tau3)
%BIAS
%    B_ = BIAS(F,IXX,IXY,IXZ,IYY,IYZ,IZZ,U1,U2,U3,U4,A,L,M,M_SAT,Q4,TAU1,TAU2,TAU3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    07-Oct-2020 07:02:58

t2 = U2.^2;
t3 = U3.^2;
t4 = m+m_sat;
t5 = 1.0./t4;
t6 = U1.^2;
t7 = a.^2;
t8 = l.^2;
t9 = q4.^2;
B_ = [t5.*(m.*tau1.*4.0+m_sat.*tau1.*4.0+Iyz.*m.*t2.*4.0-Iyz.*m.*t3.*4.0+Iyz.*m_sat.*t2.*4.0-Iyz.*m_sat.*t3.*4.0-Ixy.*U1.*U3.*m.*4.0+Ixz.*U1.*U2.*m.*4.0+Iyy.*U2.*U3.*m.*4.0-Izz.*U2.*U3.*m.*4.0-Ixy.*U1.*U3.*m_sat.*4.0+Ixz.*U1.*U2.*m_sat.*4.0+Iyy.*U2.*U3.*m_sat.*4.0-Izz.*U2.*U3.*m_sat.*4.0-U1.*U4.*a.*m.*m_sat.*8.0-U1.*U4.*l.*m.*m_sat.*4.0-U1.*U4.*m.*m_sat.*q4.*8.0+U2.*U3.*m.*m_sat.*t7.*4.0+U2.*U3.*m.*m_sat.*t8+U2.*U3.*m.*m_sat.*t9.*4.0+U2.*U3.*a.*l.*m.*m_sat.*4.0+U2.*U3.*a.*m.*m_sat.*q4.*8.0+U2.*U3.*l.*m.*m_sat.*q4.*4.0).*(-1.0./4.0);t5.*(m.*tau2.*-4.0-m_sat.*tau2.*4.0-Ixz.*m.*t3.*4.0+Ixz.*m.*t6.*4.0-Ixz.*m_sat.*t3.*4.0+Ixz.*m_sat.*t6.*4.0+Ixx.*U1.*U3.*m.*4.0-Ixy.*U2.*U3.*m.*4.0+Iyz.*U1.*U2.*m.*4.0-Izz.*U1.*U3.*m.*4.0+Ixx.*U1.*U3.*m_sat.*4.0-Ixy.*U2.*U3.*m_sat.*4.0+Iyz.*U1.*U2.*m_sat.*4.0-Izz.*U1.*U3.*m_sat.*4.0+U2.*U4.*a.*m.*m_sat.*8.0+U2.*U4.*l.*m.*m_sat.*4.0+U2.*U4.*m.*m_sat.*q4.*8.0+U1.*U3.*m.*m_sat.*t7.*4.0+U1.*U3.*m.*m_sat.*t8+U1.*U3.*m.*m_sat.*t9.*4.0+U1.*U3.*a.*l.*m.*m_sat.*4.0+U1.*U3.*a.*m.*m_sat.*q4.*8.0+U1.*U3.*l.*m.*m_sat.*q4.*4.0).*(1.0./4.0);-tau3+Ixy.*t2-Ixy.*t6-Ixx.*U1.*U2+Ixz.*U2.*U3+Iyy.*U1.*U2-Iyz.*U1.*U3;t5.*(F.*m.*2.0+F.*m_sat.*2.0+a.*m.*m_sat.*t2.*2.0+a.*m.*m_sat.*t6.*2.0+l.*m.*m_sat.*t2+l.*m.*m_sat.*t6+m.*m_sat.*q4.*t2.*2.0+m.*m_sat.*q4.*t6.*2.0).*(-1.0./2.0)];
