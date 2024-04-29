clc;
clear all %#ok<*CLALL>
simulation_time = [0:0.1:50]; %#ok<*NBRAK>
Anim = 0; % set 1 to show animation
%% Problem initializing
global m_s m_p a b c Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 tau_p1 tau_p2 px py pz
m_s = 100; Ixx = 67; Iyy = 67; Izz = 67; Ixy = 5; Ixz = 2; Iyz = 0; a = 2; b = 2; c = 2;
m_p = 10;
tau1 = 0; tau2 = 0; tau3 = 0;  tau_p1 = 0; tau_p2 = 0; 
%% Solving the system using ODE45
options = odeset('maxstep',0.1);
initial_condition = [pi/10;pi/6;0.03;  0;0;  3;3;9;      0.3;0;-0.3;  0;0.2;  0.1;-0.3;-0.4]; % Main initial condition
PXYZ = PXYZ_Voltra(b,0.3,0,-0.3,0,0.2,0.1,-0.3,-0.4,m_p,m_s,pi/10,pi/6,0.03,0,0);
uNI = U_NI(0.3,0,-0.3,0,0.2,pi/6,0.03);
px = PXYZ(1); py = PXYZ(2); pz = PXYZ(3);
z0 = [pi/10;pi/6;0.03;  0;0;  3;3;9;     uNI];
tic;
[t,z] = ode45(@Sat_Dynamics_Voltra_IC, simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7); q8 = z(:,8);
U1 = z(:,9); U2 = z(:,10); U3 = z(:,11); U4 = z(:,12); U5 = z(:,13);

%% Energy conservation criteria
E = Energy(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,U1,U2,U3,U4,U5,a,b,c,m_p,m_s,px,py,pz,q4,q5);
E_Err_Vol_IC = (E-E(1))/E(1);
% save('E_Err_Vol_IC');
x = ['Norm of Energy Error = ',num2str(norm((E-E(1))/E(1)))];
disp(x);
%% Linear momentum criteria
LM = Lmoment(px,py,pz);
LMX = LM(1)*ones(length(t));
LM_Err_Vol_IC = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Vol_IC');
x = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(x);

%% Animation 
if Anim == 1
figure;
hold on
Len=length(t);
Xcm=t; Ycm=t; Zcm=t;
psi = q1; theta = q2; phi = q3; gamma1 = q4; gamma2 = q5; X = q6; Y = q7; Z = q8;
q1 = X; q2 = Y; q3 = Z; q4 = psi; q5 = theta; q6 = phi; q7 = gamma1; q8 = gamma2;
myVideo = VideoWriter('CaseStudy2'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)
for i=1:Len
    R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];
    
    Rt=R';
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;c];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;c];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;c];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;c];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;-c];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;-c];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;-c];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;-c];
    
    
    Pg1 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(q7(i));-sin(q7(i));0]);
    Pg2 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(q8(i));sin(q8(i));0]);
      
    
    P9 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P10 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P11 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;c]+b*[-cos(q8(i));+sin(q8(i));0]);
    P12 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(q8(i));+sin(q8(i));0]);
    
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
    hold on
    box on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k','linewidth',2)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k','linewidth',2)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'b-','linewidth',1)
    plot3(V106(1,:),V106(2,:),V106(3,:),'b-','linewidth',1)
    plot3(V29(1,:),V29(2,:),V29(3,:),'b-','linewidth',1)
    plot3(V26(1,:),V26(2,:),V26(3,:),'b-','linewidth',1)
    
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',1)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',1)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',1)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',1)   
    
    
    axis equal
    axis([q1(i)- 2*a  (q1(i)+ 2*a) (q2(i)- 2*a) (q2(i)+2*a)  (q3(i)- 2*a) (q3(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(q1(i),q2(i),q3(i)-5,str)
    
    hold off
    pause(0.01)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting SnapShots
figure;
i = 1;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];
    
    Rt=R';
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;c];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;c];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;c];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;c];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;-c];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;-c];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;-c];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;-c];
    
    
    Pg1 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(q7(i));-sin(q7(i));0]);
    Pg2 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(q8(i));sin(q8(i));0]);
      
    
    P9 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P10 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P11 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;c]+b*[-cos(q8(i));+sin(q8(i));0]);
    P12 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(q8(i));+sin(q8(i));0]);
    
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
    hold on
    box on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k','linewidth',2)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k','linewidth',2)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'b-','linewidth',1)
    plot3(V106(1,:),V106(2,:),V106(3,:),'b-','linewidth',1)
    plot3(V29(1,:),V29(2,:),V29(3,:),'b-','linewidth',1)
    plot3(V26(1,:),V26(2,:),V26(3,:),'b-','linewidth',1)
    
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',1)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',1)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',1)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',1)    
    
    
    axis equal
    axis([q1(i)- 2*a  (q1(i)+ 2*a) (q2(i)- 2*a) (q2(i)+2*a)  (q3(i)- 2*a) (q3(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(5.5,8.5,str)
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 201;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];
    
    Rt=R';
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;c];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;c];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;c];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;c];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;-c];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;-c];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;-c];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;-c];
    
    
    Pg1 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(q7(i));-sin(q7(i));0]);
    Pg2 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(q8(i));sin(q8(i));0]);
      
    
    P9 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P10 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P11 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;c]+b*[-cos(q8(i));+sin(q8(i));0]);
    P12 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(q8(i));+sin(q8(i));0]);
    
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
    hold on
    box on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k','linewidth',2)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k','linewidth',2)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'b-','linewidth',1)
    plot3(V106(1,:),V106(2,:),V106(3,:),'b-','linewidth',1)
    plot3(V29(1,:),V29(2,:),V29(3,:),'b-','linewidth',1)
    plot3(V26(1,:),V26(2,:),V26(3,:),'b-','linewidth',1)
    
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',1)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',1)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',1)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',1)      
    
    
    axis equal
    axis([q1(i)- 2*a  (q1(i)+ 2*a) (q2(i)- 2*a) (q2(i)+2*a)  (q3(i)- 2*a) (q3(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(-1,-8,str)
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 401;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];
    
    Rt=R';
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;c];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;c];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;c];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;c];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;-c];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;-c];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;-c];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;-c];
    
    
    Pg1 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(q7(i));-sin(q7(i));0]);
    Pg2 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(q8(i));sin(q8(i));0]);
      
    
    P9 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P10 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P11 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;c]+b*[-cos(q8(i));+sin(q8(i));0]);
    P12 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(q8(i));+sin(q8(i));0]);
    
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
    hold on
    box on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k','linewidth',2)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k','linewidth',2)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'b-','linewidth',1)
    plot3(V106(1,:),V106(2,:),V106(3,:),'b-','linewidth',1)
    plot3(V29(1,:),V29(2,:),V29(3,:),'b-','linewidth',1)
    plot3(V26(1,:),V26(2,:),V26(3,:),'b-','linewidth',1)
    
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',1)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',1)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',1)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',1)    
    
    
    axis equal
    axis([q1(i)- 2*a  (q1(i)+ 2*a) (q2(i)- 2*a) (q2(i)+2*a)  (q3(i)- 2*a) (q3(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(-7.5,-24.5,str)
    hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 501;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];
    
    Rt=R';
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;c];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;c];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;c];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;c];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;a;-c];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[b;-a;-c];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;a;-c];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-b;-a;-c];
    
    
    Pg1 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;0]+b/2*[-cos(q7(i));-sin(q7(i));0]);
    Pg2 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;0]+b/2*[-cos(q8(i));sin(q8(i));0]);
      
    
    P9 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P10 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;-a;-c]+b*[-cos(q7(i));-sin(q7(i));0]);
    P11 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;c]+b*[-cos(q8(i));+sin(q8(i));0]);
    P12 = [q1(i);q2(i);q3(i)]+Rt*(1/2*[b;a;-c]+b*[-cos(q8(i));+sin(q8(i));0]);
    
    
    V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];
    V910 = [P9,P10]; V106 = [P10,P6]; V29 = [P2,P9];
    V411= [P4,P11];   V1112= [P11,P12]; V125= [P12,P5];
    
    plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
    hold on
    box on
    plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
    plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
    plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
    plot3(V45(1,:),V45(2,:),V45(3,:),'k','linewidth',2)
    plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
    plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
    plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
    plot3(V18(1,:),V18(2,:),V18(3,:),'k','linewidth',2)
    plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
    plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
    plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)
    
    plot3(V910(1,:),V910(2,:),V910(3,:),'b-','linewidth',1)
    plot3(V106(1,:),V106(2,:),V106(3,:),'b-','linewidth',1)
    plot3(V29(1,:),V29(2,:),V29(3,:),'b-','linewidth',1)
    plot3(V26(1,:),V26(2,:),V26(3,:),'b-','linewidth',1)
    
    plot3(V411(1,:),V411(2,:),V411(3,:),'r-','linewidth',1)
    plot3(V1112(1,:),V1112(2,:),V1112(3,:),'r-','linewidth',1)
    plot3(V125(1,:),V125(2,:),V125(3,:),'r-','linewidth',1)
    plot3(V45(1,:),V45(2,:),V45(3,:),'r-','linewidth',1)     
    
    
    axis equal
    axis([q1(i)- 2*a  (q1(i)+ 2*a) (q2(i)- 2*a) (q2(i)+2*a)  (q3(i)- 2*a) (q3(i)+2*a)])
    str=['Time = ',num2str(t(i))];
    text(-11,-33,str)
    hold off
end