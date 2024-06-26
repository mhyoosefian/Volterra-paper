function E = Energy_Voltra(P,U,g,l,m1,m2,q1,q2)
%ENERGY_VOLTRA
%    E = ENERGY_VOLTRA(P,U,G,L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    03-Oct-2020 21:52:28

t2 = q1-q2;
t3 = cos(t2);
t4 = P.^2;
t5 = q1.*2.0;
t11 = q2.*2.0;
t6 = t5-t11;
t7 = cos(t6);
t8 = U.^2;
t9 = l.^2;
t10 = m2.^2;
t12 = cos(q1);
E = -g.*l.*m1+(t4.*(3.0./2.0)+t3.*t4.*2.0+t4.*t7.*(1.0./2.0)+t8.*t9.*t10.*(1.9e1./1.6e1)+m1.*m2.*t8.*t9-t7.*t8.*t9.*t10.*(1.7e1./4.8e1)-t8.*t9.*t10.*cos(t5).*(2.5e1./3.2e1)+t8.*t9.*t10.*cos(t11).*(5.0./1.6e1)-t8.*t9.*t10.*cos(q2.*-4.0+t5).*(1.0./3.2e1)-m1.*m2.*t7.*t8.*t9.*(1.0./3.0))./((m1+m2.*2.0).*(t3.*4.0+t7+3.0))-g.*l.*m2.*t12.*(1.0./2.0)-g.*l.*m2.*(t12.*2.0+cos(q2)).*(1.0./2.0);
