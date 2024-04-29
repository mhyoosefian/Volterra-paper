clc; 
close all; 
clear all; %#ok<*CLALL>
global m1 l m2 g tau P

m1 = 1; m2 = 0.5;
l = 0.2;
tau = 0;
g = 9.81;
initcond = [pi/2; pi/2; 4; 1; -1; 3];     %[q; dq]
z0 = [pi/2;pi/2;4;  -2];     %[q; uNI]
M22 = M22_Voltra(m1,m2);
P = M22*(3);
simulation_time = 50; %#ok<*NBRAK>
option = odeset('maxstep',0.01);
tic;
[t,z] = ode45(@Dyn_Voltra, [0:0.01:simulation_time],z0,option); 
toc

q1 = z(:,1); q2 = z(:,2); q3 = z(:,3);
U = z(:,4);
dq = dq_itu_Voltra(P,U,l,m1,m2,q1,q2);
dq1 = dq(1:size(t,1),1); dq2 = dq(size(t,1)+1:2*size(t,1),1); dq3 = dq(2*size(t,1)+1:3*size(t,1),1);

E = Energy_Voltra(P,U,g,l,m1,m2,q1,q2);
E_Err_Vol_IC = (E-E(1))/E(1);
% save('E_Err_Vol_IC');
X = ['Norm of Energy Error = ',num2str(norm((E-E(1))/E(1)))];
disp(X);

f = constraints;
if f == 0
f = zeros(size(t));    %The constraint function is zero
end
f_Err_Vol_IC = f;
% save('f_Err_Vol_IC');
X = ['Norm of Constraint Error = ',num2str(norm(f))];
disp(X);

% Linear Momentum Conservation
LM = Lmom_Voltra(P,U,l,m1,m2,q1,q2);
LMX = LM(1:length(t),:);
LM_Err_Vol_IC = (LMX-LMX(1))/LMX(1);
% save('LM_Err_Vol_IC');
X = ['Norm of Momentum Error = ',num2str(norm((LMX-LMX(1))/LMX(1)))];
disp(X);

%% Animation
% close all     %Uncomment to show the animation
% figure;
% box on
% xo = 0; yo = 0;
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% for i=1:10:length(t)
%    %Points of the cart
%    xA = q3(i); yA = yo;
%    xD = xA - 2.1*l; yD = yA + 2.05*l;
%    xE = xA + 2.1*l; yE = yA + 2.05*l;
%    xG = xD; yG = yD - 4.1*l;
%    xF = xE; yF = yE - 4.1*l;
%    
%    %Points of the pendulums
%    xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
%    xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
%    
%    %Points of the wheel at the point C
%    xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
%    yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
%    xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
%    yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
%    xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
%    yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
%    xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
%    yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
%    
%    %Vectors of the cart
%    vxDE = [xD; xE]; vyDE = [yD; yE];
%    vxEF = [xE; xF]; vyEF = [yE; yF];
%    vxFG = [xF; xG]; vyFG = [yF; yG];
%    vxGD = [xG; xD]; vyGD = [yG; yD];
%    
%    %Vectors of the pendulum
%    vxAB = [xA; xB]; vyAB = [yA; yB];
%    vxBC = [xB; xC]; vyBC = [yB; yC];
%    
%    %Vectors of the wheel at point c
%    vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
%    vxJK = [xJ; xK]; vyJK = [yJ; yK];
%    vxKL = [xK; xL]; vyKL = [yK; yL];
%    vxLI = [xL; xI]; vyLI = [yL; yI];
%    
%    %Ground vector
%    vxground = [xG - 0.5*l; xF + 0.5*l];
%    vyground = [yG - 0.08; yF - 0.08];
%    
%    %Plot
%    plot(vxDE, vyDE, 'k', 'linewidth', 2)
%    hold on
%    plot(vxEF, vyEF, 'k', 'linewidth', 2)
%    plot(vxFG, vyFG, 'k', 'linewidth', 2)
%    plot(vxGD, vyGD, 'k', 'linewidth', 2)
%    plot(vxAB, vyAB, 'k', 'linewidth', 2)
%    plot(vxBC, vyBC, 'k', 'linewidth', 2)
%    plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
%    plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
%    plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
%    plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
%    plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
%    plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
%    plot(vxground , vyground, 'k', 'linewidth', 3)
%    plot(xA, yA, 'ok', 'linewidth', 5)
%    
%    str=['Time = ',num2str(t(i))];
%    text(xD + 0.3*l, yD - 0.04, str);
%    axis equal
%    axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
%    hold off
%    pause(0.01)
%    
%    frame = getframe(gcf); %get frame
%    writeVideo(myVideo, frame);
% end
% close(myVideo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting snapshots
xo = 0; yo = 0;
close all
figure;
box on
i = 1;
   %Points of the cart
   xA = q3(i); yA = yo;
   xD = xA - 2.1*l; yD = yA + 2.05*l;
   xE = xA + 2.1*l; yE = yA + 2.05*l;
   xG = xD; yG = yD - 4.1*l;
   xF = xE; yF = yE - 4.1*l;
   
   %Points of the pendulums
   xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
   xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
   
   %Points of the wheel at the point C
   xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   
   %Vectors of the cart
   vxDE = [xD; xE]; vyDE = [yD; yE];
   vxEF = [xE; xF]; vyEF = [yE; yF];
   vxFG = [xF; xG]; vyFG = [yF; yG];
   vxGD = [xG; xD]; vyGD = [yG; yD];
   
   %Vectors of the pendulum
   vxAB = [xA; xB]; vyAB = [yA; yB];
   vxBC = [xB; xC]; vyBC = [yB; yC];
   
   %Vectors of the wheel at point c
   vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
   vxJK = [xJ; xK]; vyJK = [yJ; yK];
   vxKL = [xK; xL]; vyKL = [yK; yL];
   vxLI = [xL; xI]; vyLI = [yL; yI];
   
   %Ground vector
   vxground = [xG - 0.5*l; xF + 0.5*l];
   vyground = [yG - 0.08; yF - 0.08];
   
   %Plot
   plot(vxDE, vyDE, 'k', 'linewidth', 2)
   hold on
   box on
   plot(vxEF, vyEF, 'k', 'linewidth', 2)
   plot(vxFG, vyFG, 'k', 'linewidth', 2)
   plot(vxGD, vyGD, 'k', 'linewidth', 2)
   plot(vxAB, vyAB, 'k', 'linewidth', 2)
   plot(vxBC, vyBC, 'k', 'linewidth', 2)
   plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(vxground , vyground, 'k', 'linewidth', 3)
   plot(xA, yA, 'ok', 'linewidth', 5)
   xlabel('X(m)','fontsize',10);
   ylabel('Y(m)','fontsize',10);
   str=['Time = ',num2str(t(i))];
   text(xD + 0.3*l, yD - 0.04, str);
   axis equal
   axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
   hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 1501;
   %Points of the cart
   xA = q3(i); yA = yo;
   xD = xA - 2.1*l; yD = yA + 2.05*l;
   xE = xA + 2.1*l; yE = yA + 2.05*l;
   xG = xD; yG = yD - 4.1*l;
   xF = xE; yF = yE - 4.1*l;
   
   %Points of the pendulums
   xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
   xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
   
   %Points of the wheel at the point C
   xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   
   %Vectors of the cart
   vxDE = [xD; xE]; vyDE = [yD; yE];
   vxEF = [xE; xF]; vyEF = [yE; yF];
   vxFG = [xF; xG]; vyFG = [yF; yG];
   vxGD = [xG; xD]; vyGD = [yG; yD];
   
   %Vectors of the pendulum
   vxAB = [xA; xB]; vyAB = [yA; yB];
   vxBC = [xB; xC]; vyBC = [yB; yC];
   
   %Vectors of the wheel at point c
   vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
   vxJK = [xJ; xK]; vyJK = [yJ; yK];
   vxKL = [xK; xL]; vyKL = [yK; yL];
   vxLI = [xL; xI]; vyLI = [yL; yI];
   
   %Ground vector
   vxground = [xG - 0.5*l; xF + 0.5*l];
   vyground = [yG - 0.08; yF - 0.08];
   
   %Plot
   plot(vxDE, vyDE, 'k', 'linewidth', 2)
   hold on
   box on
   plot(vxEF, vyEF, 'k', 'linewidth', 2)
   plot(vxFG, vyFG, 'k', 'linewidth', 2)
   plot(vxGD, vyGD, 'k', 'linewidth', 2)
   plot(vxAB, vyAB, 'k', 'linewidth', 2)
   plot(vxBC, vyBC, 'k', 'linewidth', 2)
   plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(vxground , vyground, 'k', 'linewidth', 3)
   plot(xA, yA, 'ok', 'linewidth', 5)
   xlabel('X(m)','fontsize',10);
   ylabel('Y(m)','fontsize',10);
   str=['Time = ',num2str(t(i))];
   text(xD + 0.3*l, yD - 0.04, str);
   axis equal
   axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
   hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 2001;
   %Points of the cart
   xA = q3(i); yA = yo;
   xD = xA - 2.1*l; yD = yA + 2.05*l;
   xE = xA + 2.1*l; yE = yA + 2.05*l;
   xG = xD; yG = yD - 4.1*l;
   xF = xE; yF = yE - 4.1*l;
   
   %Points of the pendulums
   xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
   xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
   
   %Points of the wheel at the point C
   xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   
   %Vectors of the cart
   vxDE = [xD; xE]; vyDE = [yD; yE];
   vxEF = [xE; xF]; vyEF = [yE; yF];
   vxFG = [xF; xG]; vyFG = [yF; yG];
   vxGD = [xG; xD]; vyGD = [yG; yD];
   
   %Vectors of the pendulum
   vxAB = [xA; xB]; vyAB = [yA; yB];
   vxBC = [xB; xC]; vyBC = [yB; yC];
   
   %Vectors of the wheel at point c
   vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
   vxJK = [xJ; xK]; vyJK = [yJ; yK];
   vxKL = [xK; xL]; vyKL = [yK; yL];
   vxLI = [xL; xI]; vyLI = [yL; yI];
   
   %Ground vector
   vxground = [xG - 0.5*l; xF + 0.5*l];
   vyground = [yG - 0.08; yF - 0.08];
   
   %Plot
   plot(vxDE, vyDE, 'k', 'linewidth', 2)
   hold on
   box on
   plot(vxEF, vyEF, 'k', 'linewidth', 2)
   plot(vxFG, vyFG, 'k', 'linewidth', 2)
   plot(vxGD, vyGD, 'k', 'linewidth', 2)
   plot(vxAB, vyAB, 'k', 'linewidth', 2)
   plot(vxBC, vyBC, 'k', 'linewidth', 2)
   plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(vxground , vyground, 'k', 'linewidth', 3)
   plot(xA, yA, 'ok', 'linewidth', 5)
   xlabel('X(m)','fontsize',10);
   ylabel('Y(m)','fontsize',10);
   str=['Time = ',num2str(t(i))];
   text(xD + 0.3*l, yD - 0.04, str);
   axis equal
   axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
   hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 3501;
   %Points of the cart
   xA = q3(i); yA = yo;
   xD = xA - 2.1*l; yD = yA + 2.05*l;
   xE = xA + 2.1*l; yE = yA + 2.05*l;
   xG = xD; yG = yD - 4.1*l;
   xF = xE; yF = yE - 4.1*l;
   
   %Points of the pendulums
   xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
   xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
   
   %Points of the wheel at the point C
   xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   
   %Vectors of the cart
   vxDE = [xD; xE]; vyDE = [yD; yE];
   vxEF = [xE; xF]; vyEF = [yE; yF];
   vxFG = [xF; xG]; vyFG = [yF; yG];
   vxGD = [xG; xD]; vyGD = [yG; yD];
   
   %Vectors of the pendulum
   vxAB = [xA; xB]; vyAB = [yA; yB];
   vxBC = [xB; xC]; vyBC = [yB; yC];
   
   %Vectors of the wheel at point c
   vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
   vxJK = [xJ; xK]; vyJK = [yJ; yK];
   vxKL = [xK; xL]; vyKL = [yK; yL];
   vxLI = [xL; xI]; vyLI = [yL; yI];
   
   %Ground vector
   vxground = [xG - 0.5*l; xF + 0.5*l];
   vyground = [yG - 0.08; yF - 0.08];
   
   %Plot
   plot(vxDE, vyDE, 'k', 'linewidth', 2)
   hold on
   box on
   plot(vxEF, vyEF, 'k', 'linewidth', 2)
   plot(vxFG, vyFG, 'k', 'linewidth', 2)
   plot(vxGD, vyGD, 'k', 'linewidth', 2)
   plot(vxAB, vyAB, 'k', 'linewidth', 2)
   plot(vxBC, vyBC, 'k', 'linewidth', 2)
   plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(vxground , vyground, 'k', 'linewidth', 3)
   plot(xA, yA, 'ok', 'linewidth', 5)
   xlabel('X(m)','fontsize',10);
   ylabel('Y(m)','fontsize',10);
   str=['Time = ',num2str(t(i))];
   text(xD + 0.3*l, yD - 0.04, str);
   axis equal
   axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
   hold off
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
i = 5001;
   %Points of the cart
   xA = q3(i); yA = yo;
   xD = xA - 2.1*l; yD = yA + 2.05*l;
   xE = xA + 2.1*l; yE = yA + 2.05*l;
   xG = xD; yG = yD - 4.1*l;
   xF = xE; yF = yE - 4.1*l;
   
   %Points of the pendulums
   xB = xA + l*sin(q1(i)); yB = yA - l*cos(q1(i));
   xC = xB + l*sin(q2(i)); yC = yB - l*cos(q2(i));
   
   %Points of the wheel at the point C
   xI = xB + 0.9*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yI = yB - 0.9*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xJ = xB + 1.1*l*sin(q2(i)) + 0.05*l*cos(q2(i)); 
   yJ = yB - 1.1*l*cos(q2(i)) + 0.05*l*sin(q2(i));
   xL = xB + 0.9*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yL = yB - 0.9*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   xK = xB + 1.1*l*sin(q2(i)) - 0.05*l*cos(q2(i)); 
   yK = yB - 1.1*l*cos(q2(i)) - 0.05*l*sin(q2(i));
   
   %Vectors of the cart
   vxDE = [xD; xE]; vyDE = [yD; yE];
   vxEF = [xE; xF]; vyEF = [yE; yF];
   vxFG = [xF; xG]; vyFG = [yF; yG];
   vxGD = [xG; xD]; vyGD = [yG; yD];
   
   %Vectors of the pendulum
   vxAB = [xA; xB]; vyAB = [yA; yB];
   vxBC = [xB; xC]; vyBC = [yB; yC];
   
   %Vectors of the wheel at point c
   vxIJ = [xI; xJ]; vyIJ = [yI; yJ];
   vxJK = [xJ; xK]; vyJK = [yJ; yK];
   vxKL = [xK; xL]; vyKL = [yK; yL];
   vxLI = [xL; xI]; vyLI = [yL; yI];
   
   %Ground vector
   vxground = [xG - 0.5*l; xF + 0.5*l];
   vyground = [yG - 0.08; yF - 0.08];
   
   %Plot
   plot(vxDE, vyDE, 'k', 'linewidth', 2)
   hold on
   box on
   plot(vxEF, vyEF, 'k', 'linewidth', 2)
   plot(vxFG, vyFG, 'k', 'linewidth', 2)
   plot(vxGD, vyGD, 'k', 'linewidth', 2)
   plot(vxAB, vyAB, 'k', 'linewidth', 2)
   plot(vxBC, vyBC, 'k', 'linewidth', 2)
   plot(vxIJ, vyIJ, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxJK, vyJK, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxKL, vyKL, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(vxLI, vyLI, 'k', 'linewidth', 2,'color', [0.5 0.5 0.5])
   plot(xG + 0.5*l, yG - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(xF - 0.5*l, yF - 0.04, 'o','color', [0.5 0.5 0.5], 'linewidth', 18)
   plot(vxground , vyground, 'k', 'linewidth', 3)
   plot(xA, yA, 'ok', 'linewidth', 5)
   xlabel('X(m)','fontsize',10);
   ylabel('Y(m)','fontsize',10);
   str=['Time = ',num2str(t(i))];
   text(xD + 0.3*l, yD - 0.04, str);
   axis equal
   axis ([xD-0.5*l xE+0.5*l yG-0.4*l yD+0.1*l])
   hold off