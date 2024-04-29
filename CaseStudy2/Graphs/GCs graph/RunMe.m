
% Firstly, manually load all .mat files in the folder and then run the code

figure; 
hold on;
box on
plot(t, q1_Lag, 'r', 'linewidth', 2);
plot(t, q2_Lag, 'b', 'linewidth', 2);
plot(t, q3_Lag, 'g', 'linewidth', 2);
plot(t, q4_Lag, 'k', 'linewidth', 2);
plot(t, q5_Lag, 'y', 'linewidth', 2);
plot(t, q6_Lag, 'r--', 'linewidth', 2);
plot(t, q7_Lag, 'b--', 'linewidth', 2);
plot(t, q8_Lag, 'k--', 'linewidth', 2);
legend('\psi','\theta','\phi','\gamma_1','\gamma_2','X','Y', 'Z', 'orientation', 'horizontal');
xlabel('Time(s)','fontsize',12);
grid on;
set(gca,'fontsize',12);









% figure;
% hold on
% box on;
% plot(t,q1_Lag,'k','linewidth',2);
% plot(t,q1_Gibbs,'ob:','linewidth',2);
% plot(t,q1_Mag,'*y--','linewidth',2);
% plot(t,q1_Vol_IC,'r--','linewidth',2);
% legend('Lagrange','Gibbs-Appell(Kane)','Maggi','Proposed');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('\psi','fontsize',10);


% figure;
% hold on
% box on;
% plot(t,q3_Lag,'k','linewidth',2);
% plot(t,q3_Gibbs,'ob:','linewidth',2);
% plot(t,q3_Mag,'*y--','linewidth',2);
% plot(t,q3_Vol_IC,'r--','linewidth',2);
% legend('Lagrange','Gibbs-Appell(Kane)','Maggi','Proposed');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('\phi','fontsize',10);

% figure;
% hold on
% box on;
% plot(t,q5_Lag,'k','linewidth',2);
% plot(t,q5_Gibbs,'b-','linewidth',2);
% plot(t,q5_Mag,'y--','linewidth',2);
% plot(t,q5_Vol_IC,'r--','linewidth',2);
% legend('Lagrange','Gibbs-Appell(Kane)','Maggi','Proposed');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('\gamma_2','fontsize',10);

% figure;
% hold on
% box on;
% plot(t,q8_Lag,'k','linewidth',2);
% plot(t,q8_Gibbs,'ob:','linewidth',2);
% plot(t,q8_Mag,'*y--','linewidth',2);
% plot(t,q8_Vol_IC,'r--','linewidth',2);
% legend('Lagrange','Gibbs-Appell(Kane)','Maggi','Proposed');
% set(gca,'fontsize',10,'fontweight','bold');
% xlabel('Time (s)','fontsize',10);
% ylabel('Z','fontsize',10);