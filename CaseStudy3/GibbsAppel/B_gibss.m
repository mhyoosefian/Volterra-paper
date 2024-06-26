function B_g = B_gibss(F,Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,a,l,m,q1,q2,q3,q4,tau1,tau2,tau3)
%B_GIBSS
%    B_G = B_GIBSS(F,IXX,IXY,IXZ,IYY,IYZ,IZZ,U1,U2,U3,U4,A,L,M,Q1,Q2,Q3,Q4,TAU1,TAU2,TAU3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    07-Oct-2020 06:53:25

t2 = U3.^2;
t3 = a.^2;
t4 = l.^2;
t5 = q4.^2;
t6 = U1.^2;
t7 = U2.^2;
t8 = sin(q1);
t9 = sin(q3);
t10 = cos(q1);
t11 = cos(q3);
t12 = sin(q2);
t13 = cos(q2);
B_g = [-tau1+Iyz.*t2-Iyz.*t7+Ixy.*U1.*U3-Ixz.*U1.*U2-Iyy.*U2.*U3+Izz.*U2.*U3+U1.*U4.*a.*m.*2.0+U1.*U4.*l.*m+U1.*U4.*m.*q4.*2.0-U2.*U3.*m.*t3-U2.*U3.*m.*t4.*(1.0./4.0)-U2.*U3.*m.*t5-U2.*U3.*a.*l.*m-U2.*U3.*a.*m.*q4.*2.0-U2.*U3.*l.*m.*q4;-tau2-Ixz.*t2+Ixz.*t6+Ixx.*U1.*U3-Ixy.*U2.*U3+Iyz.*U1.*U2-Izz.*U1.*U3+U2.*U4.*a.*m.*2.0+U2.*U4.*l.*m+U2.*U4.*m.*q4.*2.0+U1.*U3.*m.*t3+U1.*U3.*m.*t4.*(1.0./4.0)+U1.*U3.*m.*t5+U1.*U3.*a.*l.*m+U1.*U3.*a.*m.*q4.*2.0+U1.*U3.*l.*m.*q4;-tau3-Ixy.*t6+Ixy.*t7-Ixx.*U1.*U2+Ixz.*U2.*U3+Iyy.*U1.*U2-Iyz.*U1.*U3;-F-a.*m.*t6-a.*m.*t7-l.*m.*t6.*(1.0./2.0)-l.*m.*t7.*(1.0./2.0)-m.*q4.*t6-m.*q4.*t7;U1.*U4.*m.*t8.*t11.*2.0+U2.*U4.*m.*t10.*t13.*2.0-a.*m.*t6.*t8.*t9-a.*m.*t7.*t8.*t9-l.*m.*t6.*t8.*t9.*(1.0./2.0)-l.*m.*t7.*t8.*t9.*(1.0./2.0)-m.*q4.*t6.*t8.*t9-m.*q4.*t7.*t8.*t9-l.*m.*t6.*t10.*t11.*t12.*(1.0./2.0)-l.*m.*t7.*t10.*t11.*t12.*(1.0./2.0)-m.*q4.*t6.*t10.*t11.*t12-m.*q4.*t7.*t10.*t11.*t12-U2.*U3.*a.*m.*t8.*t11+U1.*U3.*a.*m.*t10.*t13-U2.*U3.*l.*m.*t8.*t11.*(1.0./2.0)+U1.*U3.*l.*m.*t10.*t13.*(1.0./2.0)-U2.*U3.*m.*q4.*t8.*t11+U1.*U3.*m.*q4.*t10.*t13-U1.*U4.*m.*t9.*t10.*t12.*2.0-a.*m.*t6.*t10.*t11.*t12-a.*m.*t7.*t10.*t11.*t12+U2.*U3.*a.*m.*t9.*t10.*t12+U2.*U3.*l.*m.*t9.*t10.*t12.*(1.0./2.0)+U2.*U3.*m.*q4.*t9.*t10.*t12;U1.*U4.*m.*t10.*t11.*-2.0+U2.*U4.*m.*t8.*t13.*2.0+a.*m.*t6.*t9.*t10+a.*m.*t7.*t9.*t10+l.*m.*t6.*t9.*t10.*(1.0./2.0)+l.*m.*t7.*t9.*t10.*(1.0./2.0)+m.*q4.*t6.*t9.*t10+m.*q4.*t7.*t9.*t10-l.*m.*t6.*t8.*t11.*t12.*(1.0./2.0)-l.*m.*t7.*t8.*t11.*t12.*(1.0./2.0)-m.*q4.*t6.*t8.*t11.*t12-m.*q4.*t7.*t8.*t11.*t12+U1.*U3.*a.*m.*t8.*t13+U2.*U3.*a.*m.*t10.*t11+U1.*U3.*l.*m.*t8.*t13.*(1.0./2.0)+U2.*U3.*l.*m.*t10.*t11.*(1.0./2.0)+U1.*U3.*m.*q4.*t8.*t13+U2.*U3.*m.*q4.*t10.*t11-U1.*U4.*m.*t8.*t9.*t12.*2.0-a.*m.*t6.*t8.*t11.*t12-a.*m.*t7.*t8.*t11.*t12+U2.*U3.*a.*m.*t8.*t9.*t12+U2.*U3.*l.*m.*t8.*t9.*t12.*(1.0./2.0)+U2.*U3.*m.*q4.*t8.*t9.*t12;m.*(U2.*U4.*t12.*4.0+a.*t6.*t11.*t13.*2.0+a.*t7.*t11.*t13.*2.0+l.*t6.*t11.*t13+l.*t7.*t11.*t13+q4.*t6.*t11.*t13.*2.0+q4.*t7.*t11.*t13.*2.0+U1.*U3.*a.*t12.*2.0+U1.*U3.*l.*t12+U1.*U3.*q4.*t12.*2.0+U1.*U4.*t9.*t13.*4.0-U2.*U3.*a.*t9.*t13.*2.0-U2.*U3.*l.*t9.*t13-U2.*U3.*q4.*t9.*t13.*2.0).*(-1.0./2.0)];
