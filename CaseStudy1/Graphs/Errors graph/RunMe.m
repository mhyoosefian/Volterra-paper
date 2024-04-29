
% Firstly, manually load all .mat files in the folder and then run the code

figure;
hold on
box on;
plot(t,100*E_Err_Vol_IC,'r','linewidth',2);
plot(t,100*E_Err_Lag,'k--','linewidth',2);
plot(t,100*E_Err_Gibbs,'b-','linewidth',2);
plot(t,100*E_Err_Magstd,'y--','linewidth',2);
legend('Proposed','Lagrange','Gibbs-Appell (Kane)','Maggi');
set(gca,'fontsize',12);
grid on;
xlabel('Time (s)','fontsize',12);
ylabel('% Error of  Mechanical Energy','fontsize',12);


figure;
hold on
box on;
plot(t,100*f_Err_Vol_IC(1:length(t)),'r','linewidth',2);
plot(t,100*f_Err_Lag(1:length(t)),'k--','linewidth',2);
plot(t,100*f_Err_Gibbs(1:length(t)),'b--','linewidth',2);
plot(t,100*f_Err_Magstd(1:length(t)),'y--','linewidth',2);
legend('Proposed','Lagrange','Gibbs-Appell (Kane)','Maggi');
set(gca,'fontsize',12);
grid on;
xlabel('Time (s)','fontsize',12);
ylabel('% Error of  Constraint Conservation','fontsize',12);


figure;
hold on
box on;
plot(t,100*LM_Err_Vol_IC,'r','linewidth',2);
plot(t,100*LM_Err_Lag,'k--','linewidth',2);
plot(t,100*LM_Err_Gibbs,'b-','linewidth',2);
plot(t,100*LM_Err_Magstd,'y--','linewidth',2);
legend('Proposed','Lagrange','Gibbs-Appell (Kane)','Maggi');
set(gca,'fontsize',12);
grid on;
xlabel('Time (s)','fontsize',12);
ylabel('% Error of  Linear Momentum (X direction) ','fontsize',12);