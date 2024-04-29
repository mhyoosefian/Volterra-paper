
% Firstly, manually load all .mat files in the folder and then run the code

figure;
hold on
box on;
plot(t,100*Pow_Err_Lag,'k','linewidth',2);
plot(t,100*Pow_Err_Gibbs,'b-','linewidth',2);
plot(t,100*Pow_Err_Mag,'y--','linewidth',2);
plot(t,100*Pow_Err_Vol,'r--','linewidth',2);
legend('Lagrange','Gibbs-Appell (Kane)','Maggi','Proposed');
set(gca,'fontsize',12);
grid on;
xlabel('Time (s)','fontsize',12);
ylabel('% Energy Error ','fontsize',12);


figure;
hold on
box on;
plot(t,100*LM_Err_Lag,'k','linewidth',2);
plot(t,100*LM_Err_Gibbs,'b-','linewidth',2);
plot(t,100*LM_Err_Mag,'y','linewidth',2);
plot(t,100*LM_Err_Vol_IC,'r--','linewidth',2);
legend('Lagrange','Gibbs-Appell (Kane)','Maggi','Proposed');
set(gca,'fontsize',12);
grid on;
xlabel('Time (s)','fontsize',10);
ylabel('% Error of  Linear Momentum (X direction) ','fontsize',12);