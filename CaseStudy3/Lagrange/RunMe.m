clc;
close all;
clear all; %#ok<*CLALL>

%% Choosing whether to see Animation or not
Anim = 0; % Set 1 to show the animation
simulation_time = [0:0.1:50];
%% Problem initializing
global m_sat l Ixx Iyy Izz Ixy Ixz Iyz tau1 tau2 tau3 a m
m_sat = 2000; Ixx = 1400; Iyy = 900; Izz = 1100; Ixy = 5; Ixz = 8; Iyz = 3; l = 2.5;
m = 1; a = 0.5;
tau1 = 0; tau2 = 0; tau3 = 0;
pow_0 = Power(0.012,-0.1,0.05,-0.05,0,pi/12,0.08,tau1,tau2,tau3);
z0 = [pi/8;pi/12;0.08; a; 0;0;0;  -0.1;0.05;-0.05; 0; 2;1;0;    pow_0]; % Initial condition

%% Solving the system using ODE45
options = odeset('maxstep',0.1);
tic;
[t,z] = ode45(@Sat_Dynamics_Lagrange,simulation_time,z0,options);  %#ok<*NBRAK>
toc
q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6); q7 = z(:,7);
dq1 = z(:,8); dq2 = z(:,9); dq3 = z(:,10); dq4 = z(:,11); dq5 = z(:,12); dq6 = z(:,13); dq7 = z(:,14);
q1_Lag = q1; q2_Lag = q2; q3_Lag = q3; q4_Lag = q4; q5_Lag = q5; q6_Lag = q6; q7_Lag = q7;
% save ('q1_Lag'); save ('q2_Lag'); save ('q3_Lag'); save ('q4_Lag'); save ('q5_Lag'); save ('q6_Lag'); save ('q7_Lag');
int_pow = z(:,15);

%% Power conservation criteria
if Anim == 0
F = -0.018*sin(0.089*t) + 0.012*cos(0.0485*t);
E = Energy(Ixx,Ixy,Ixz,Iyy,Iyz,Izz,a,dq1,dq2,dq3,dq4,dq5,dq6,dq7,l,m,m_sat,q1,q2,q3,q4);
Pow_Err_Lag = (E-int_pow);
Pow_Err_Lag = (Pow_Err_Lag - Pow_Err_Lag(1))/Pow_Err_Lag(1);
% save('Pow_Err_Lag');
X = ['Norm of Power Error = ',num2str(norm(Pow_Err_Lag))];
disp(X);
end

%% Linear momentum criteria
if Anim == 0
    LM = Lmoment(a,dq1,dq2,dq3,dq4,dq5,dq6,dq7,l,m,m_sat,q1,q2,q3,q4);
    LMX = LM(1:length(t));
    LMY = LM(length(t)+1:2*length(t));
    LMZ = LM(2*length(t)+1:3*length(t));
    LM_Err_Lag = (LMX-LMX(1))/LMX(1);
%     save('LM_Err_Lag');
    X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
    disp(X);
end

%% Animation 
if  Anim == 1
q1 = q5_Lag; q2 = q6_Lag; q3 = q7_Lag; q4 = q1_Lag; q5 = q2_Lag; q6 = q3_Lag; q7 = q4_Lag;
figure
box on
hold on
Len=length(t);
myVideo = VideoWriter('CaseStudy3'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)
for i=1:Len
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];

Rt = R';

% Points of Satelite
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;l];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;l];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;l];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;l];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;-l];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;-l];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;-l];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;-l];

P9 = [q1(i);q2(i);q3(i)]+Rt*1/2*[0;0;l];
P10 = [q1(i);q2(i);q3(i)]+Rt*[0;0;l/2+a+q7(i)];

V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];


V910 = [P9,P10];


plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
box on
hold on
plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
plot3(V45(1,:),V45(2,:),V45(3,:),'r','linewidth',2)
plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
plot3(V18(1,:),V18(2,:),V18(3,:),'b','linewidth',2)
plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)

plot3(V910(1,:),V910(2,:),V910(3,:),'k','linewidth',1)

plot3(P10(1), P10(2),P10(3), 'ok', 'linewidth', 4)

axis equal
axis([q1(i)- 2*l  (q1(i)+ 3*l) (q2(i)- 2*l) (q2(i)+3*l)  (q3(i)- 2*l) (q3(i)+3*l)])
str=['Time = ',num2str(t(i))];
text(q1(i),q2(i),q3(i)-5,str)

hold off
pause(0.01)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Snap Shots
figure;
i = 1;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];

Rt = R';
% Points of Satelite
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;l];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;l];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;l];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;l];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;-l];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;-l];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;-l];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;-l];

P9 = [q1(i);q2(i);q3(i)]+Rt*1/2*[0;0;l];
P10 = [q1(i);q2(i);q3(i)]+Rt*[0;0;l/2+a+q7(i)];

V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];


V910 = [P9,P10];


plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
box on
hold on
plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
plot3(V45(1,:),V45(2,:),V45(3,:),'r','linewidth',2)
plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
plot3(V18(1,:),V18(2,:),V18(3,:),'b','linewidth',2)
plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)

plot3(V910(1,:),V910(2,:),V910(3,:),'k','linewidth',1)

plot3(P10(1), P10(2),P10(3), 'ok', 'linewidth', 4)

axis equal
axis([q1(i)- 2*l  (q1(i)+ 3*l) (q2(i)- 2*l) (q2(i)+3*l)  (q3(i)- 2*l) (q3(i)+3*l)])
str=['Time = ',num2str(t(i))];
text(-9,-8,str);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 176;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];

Rt = R';
% Points of Satelite
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;l];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;l];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;l];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;l];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;-l];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;-l];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;-l];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;-l];

P9 = [q1(i);q2(i);q3(i)]+Rt*1/2*[0;0;l];
P10 = [q1(i);q2(i);q3(i)]+Rt*[0;0;l/2+a+q7(i)];

V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];


V910 = [P9,P10];


plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
box on
hold on
plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
plot3(V45(1,:),V45(2,:),V45(3,:),'r','linewidth',2)
plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
plot3(V18(1,:),V18(2,:),V18(3,:),'b','linewidth',2)
plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)

plot3(V910(1,:),V910(2,:),V910(3,:),'k','linewidth',1)

plot3(P10(1), P10(2),P10(3), 'ok', 'linewidth', 4)

axis equal
axis([q1(i)- 2*l  (q1(i)+ 3*l) (q2(i)- 2*l) (q2(i)+3*l)  (q3(i)- 2*l) (q3(i)+3*l)])
str=['Time = ',num2str(t(i))];
text(26,10,str);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 351;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];

Rt = R';

% Points of Satelite
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;l];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;l];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;l];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;l];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;-l];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;-l];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;-l];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;-l];

P9 = [q1(i);q2(i);q3(i)]+Rt*1/2*[0;0;l];
P10 = [q1(i);q2(i);q3(i)]+Rt*[0;0;l/2+a+q7(i)];

V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];


V910 = [P9,P10];


plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
box on
hold on
plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
plot3(V45(1,:),V45(2,:),V45(3,:),'r','linewidth',2)
plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
plot3(V18(1,:),V18(2,:),V18(3,:),'b','linewidth',2)
plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)

plot3(V910(1,:),V910(2,:),V910(3,:),'k','linewidth',1)

plot3(P10(1), P10(2),P10(3), 'ok', 'linewidth', 4)

axis equal
axis([q1(i)- 2*l  (q1(i)+ 3*l) (q2(i)- 2*l) (q2(i)+3*l)  (q3(i)- 2*l) (q3(i)+3*l)])
str=['Time = ',num2str(t(i))];
text(61,27,str);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 501;
R = [                          cos(q4(i))*cos(q5(i)),                           cos(q5(i))*sin(q4(i)),        -sin(q5(i))
        cos(q4(i))*sin(q5(i))*sin(q6(i)) - cos(q6(i))*sin(q4(i)), cos(q4(i))*cos(q6(i)) + sin(q4(i))*sin(q5(i))*sin(q6(i)), cos(q5(i))*sin(q6(i))
        sin(q4(i))*sin(q6(i)) + cos(q4(i))*cos(q6(i))*sin(q5(i)), cos(q6(i))*sin(q4(i))*sin(q5(i)) - cos(q4(i))*sin(q6(i)), cos(q5(i))*cos(q6(i))];

Rt = R';

% Points of Satelite
    P1 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;l];
    P2 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;l];
    P3 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;l];
    P4 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;l];
    P5 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;l;-l];
    P6 = [q1(i);q2(i);q3(i)]+Rt*1/2*[l;-l;-l];
    P7 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;l;-l];
    P8 = [q1(i);q2(i);q3(i)]+Rt*1/2*[-l;-l;-l];

P9 = [q1(i);q2(i);q3(i)]+Rt*1/2*[0;0;l];
P10 = [q1(i);q2(i);q3(i)]+Rt*[0;0;l/2+a+q7(i)];

V12 = [P1,P2]; V13 = [P1,P3]; V24 = [P2,P4]; V34 = [P3,P4]; V45 = [P4,P5]; V56 = [P5,P6]; V57 = [P5,P7]; V78 = [P7,P8]; V18 = [P1,P8]; V37 = [P3,P7]; V26 = [P2,P6]; V86 = [P8,P6];


V910 = [P9,P10];


plot3(V12(1,:),V12(2,:),V12(3,:),'k','linewidth',2)
box on
hold on
plot3(V13(1,:),V13(2,:),V13(3,:),'k','linewidth',2)
plot3(V24(1,:),V24(2,:),V24(3,:),'k','linewidth',2)
plot3(V34(1,:),V34(2,:),V34(3,:),'k','linewidth',2)
plot3(V45(1,:),V45(2,:),V45(3,:),'r','linewidth',2)
plot3(V56(1,:),V56(2,:),V56(3,:),'k','linewidth',2)
plot3(V57(1,:),V57(2,:),V57(3,:),'k','linewidth',2)
plot3(V78(1,:),V78(2,:),V78(3,:),'k','linewidth',2)
plot3(V18(1,:),V18(2,:),V18(3,:),'b','linewidth',2)
plot3(V37(1,:),V37(2,:),V37(3,:),'k','linewidth',2)
plot3(V26(1,:),V26(2,:),V26(3,:),'k','linewidth',2)
plot3(V86(1,:),V86(2,:),V86(3,:),'k','linewidth',2)

plot3(V910(1,:),V910(2,:),V910(3,:),'k','linewidth',1)

plot3(P10(1), P10(2),P10(3), 'ok', 'linewidth', 4)

axis equal
axis([q1(i)- 2*l  (q1(i)+ 3*l) (q2(i)- 2*l) (q2(i)+3*l)  (q3(i)- 2*l) (q3(i)+3*l)])
str=['Time = ',num2str(t(i))];
text(91,42,str);
hold off
end