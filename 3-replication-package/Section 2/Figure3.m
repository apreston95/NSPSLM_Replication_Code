clear all
close all
clc

%% Parameters

mu = 6;
rho = 0.05;
chi = 0.7;
gamma = 10;
wbar = 0.62;
beta = 0.965;
kappa = 1;


theta_0 = -wbar/kappa;
theta_1 = (1-(1/mu))/kappa;

lambda_0 = beta*(1 + rho*(chi^(-gamma) - 1)*(1-theta_0));
lambda_1 = beta*rho*(chi^(-gamma) - 1)*theta_1;

%% Variables


A_1 = 0.9:0.005:1.1;

Eta_1 = theta_0 + theta_1.*A_1;

M_1 = lambda_0 - lambda_1.*A_1;

R_0 = 1./M_1;

Eta_0 = (1./kappa).*(  ((mu-1)/mu) - wbar + (1-rho)*kappa.*M_1.*Eta_1);

N_0 = 1 - rho + rho.*Eta_0;

Y_0 = N_0;


V_0 = rho.*Eta_0.^2;

V_1 = rho.*Eta_1.^2;

e_1 = 1 - N_0 + rho*N_0;

N_1 = (1-rho).* N_0 + Eta_1.*e_1;

K = 0.5*(length(A_1) - 1) + 1;

Mean_Eta1 = Eta_1(K)

Mean_Eta0 = Eta_0(K);

Mean_Y0 = Y_0(K);

Mean_R0 = R_0(K);

Mean_N1 = N_1(K);









%% Repeat when shutting off PS channel

chi=1;

theta_0 = -wbar/kappa;
theta_1 = (1-(1/mu))/kappa;

lambda_0 = beta*(1 + rho*(chi^(-gamma) - 1)*(1-theta_0));
lambda_1 = beta*rho*(chi^(-gamma) - 1)*theta_1;

Eta_1_Star = theta_0 + theta_1.*A_1;

M_1_Star = lambda_0 - lambda_1.*A_1;

R_0_Star = 1./M_1_Star;

Eta_0_Star = (1./kappa).*(  ((mu-1)/mu) - wbar + (1-rho)*beta*kappa.*M_1_Star.*Eta_1_Star);

N_0_Star = 1 - rho + rho.*Eta_0_Star;

Y_0_Star = N_0_Star;

V_0_Star = rho.*Eta_0_Star.^2;




%% Plot of Y_0 vs A_1

figure;
%plot(log(A_1), log(Y_0/Mean_Y0),'color','black','LineWidth',2);
plot(100.*log(A_1), 100.*log(Y_0./Y_0(K)),'color','red','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(Y_0_Star./Y_0_Star(K)),'color','blue','LineWidth',2,'LineStyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{Y}_0$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
title('Output', 'interpreter', 'latex', 'FontSize', 16);

grid on

print('Figure3','-dpng');