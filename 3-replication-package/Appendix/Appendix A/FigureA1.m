clear all
clc

A_1_Grid = 0.9:0.005:1.1;

for i = 1:length(A_1_Grid)

par.A_1 = A_1_Grid(i);

par.chi = 0.79;

save params par;


dynare Two_Period_Dynare_BondsOnly nolog;

Y_0_Grid(i) = oo_.steady_state(20);

B_0_E_Grid(i) = oo_.steady_state(14);

end

for i = 1:length(A_1_Grid)

par.A_1 = A_1_Grid(i);

par.chi = 1;

save params par;


dynare Two_Period_Dynare_BondsOnly nolog;


Y_0_Star(i) = oo_.steady_state(20);

B_0_E_Star(i) = oo_.steady_state(14);

end


%% 
close all

load('Y_0.mat')
K = 21;

figure;
plot(100.*log(A_1_Grid), 100.*log(Y_0./Y_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1_Grid), 100.*log(Y_0_Grid./Y_0_Grid(K)),'color','red','LineWidth',2,'linestyle','--');
hold on
plot(100.*log(A_1_Grid), 100*log(Y_0_Star./Y_0_Star(K)),'color','blue','LineWidth',2,'linestyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{Y}_0$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
grid on

print('FigureA1','-dpng');

% figure;
% plot(100.*log(A_1_Grid), 100.*log(B_0_E_Grid./B_0_E_Grid(K)),'color','black','LineWidth',2);
% hold on
% xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
% ylabel('$\hat{B}_{0}^{E}$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
% xlim( [-10 10])
% grid on
% 
% print('B_0_E_Plot_NonDegen','-dpng');