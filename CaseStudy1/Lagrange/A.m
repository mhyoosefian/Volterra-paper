function A_IM = A(l,m1,m2,q1,q2)
%A
%    A_IM = A(L,M1,M2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    02-Oct-2020 17:51:12

t2 = l.^2;
t3 = q1-q2;
t4 = cos(t3);
t5 = m2.*t2.*t4.*(1.0./2.0);
t6 = cos(q1);
t7 = l.*m2.*t6.*(3.0./2.0);
t8 = cos(q2);
t9 = l.*m2.*t8.*(1.0./2.0);
A_IM = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,m2.*t2.*(4.0./3.0),t5,t7,-l.*t4,0.0,0.0,0.0,t5,m2.*t2.*(1.0./3.0),t9,-l,0.0,0.0,0.0,t7,t9,m1+m2.*2.0,0.0,0.0,0.0,0.0,-l.*t4,-l,0.0,0.0],[7,7]);
