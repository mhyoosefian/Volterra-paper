
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
legend('\psi','\theta','\phi','\rho','X','Y', 'Z', 'orientation', 'horizontal');
xlabel('Time(s)','fontsize',12);
set(gca,'fontsize',12);
grid on;