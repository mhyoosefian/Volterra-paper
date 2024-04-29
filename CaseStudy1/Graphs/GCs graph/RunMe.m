
% Firstly, manually load all .mat files in the folder and then run the code

figure
hold on;
box on;
plot(t,q1_Lag,'r-','linewidth',2)
plot(t,q2_Lag,'b-','linewidth',2)
plot(t,q3_Lag,'g-','linewidth',2)
legend('\theta_1 (rad)','\theta_2 (rad)', 'x (m)')
xlabel('Time(s)','fontsize',12);
grid on;
set(gca,'fontsize',12);














% figure;
% hold on
% box on;
% plot(t,q1_Vol_IC,'r','linewidth',2);
% plot(t,q1_Lag,'k:','linewidth',2);
% plot(t,q1_Gibbs,'ob','linewidth',2);
% plot(t,q1_Mag,'*y','linewidth',2);
% legend('Proposed','Lagrange','Gibbs-Appell(Kane)','Maggi');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('\theta_1','fontsize',10);
% 
% 
% figure;
% hold on
% box on;
% plot(t,q2_Vol_IC(1:length(t)),'r','linewidth',2);
% plot(t,q2_Lag(1:length(t)),'k:','linewidth',2);
% plot(t,q2_Gibbs(1:length(t)),'ob','linewidth',2);
% plot(t,q2_Mag(1:length(t)),'*y','linewidth',2);
% legend('Proposed','Lagrange','Gibbs-Appell(Kane)','Maggi');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('\theta_2','fontsize',10);
% 
% 
% figure;
% hold on
% box on;
% plot(t,q3_Vol_IC,'r','linewidth',2);
% plot(t,q3_Lag,'k:','linewidth',2);
% plot(t,q3_Gibbs,'ob','linewidth',2);
% plot(t,q3_Mag,'*y','linewidth',2);
% legend('Proposed','Lagrange','Gibbs-Appell(Kane)','Maggi');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('X','fontsize',10);