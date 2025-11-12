clear all
close all
clc

%% Parameters

kappa = 1;
psi = 200;
mu = 6;
rho = 0.05;
chi = 0.79;
gamma = 3;
wbar = 0.52;
beta = 0.965;

delta = 1;

theta_0 = -wbar/kappa;
theta_1 = (1-(1/mu))/kappa;

lambda_0 = beta*(1 + rho*(chi^(-gamma) - 1)*(1-theta_0));
lambda_1 = beta*rho*(chi^(-gamma) - 1)*theta_1;

Rbar = 1;

%% Variables


A_1 = 0.9:0.005:1.1;

Eta_1 = theta_0 + theta_1.*A_1;

M_1 = lambda_0 - lambda_1.*A_1;

R_0 = 1./M_1;

P_0 = (R_0/Rbar).^(1/delta);

Eta_0 = (1./kappa).*( (psi./mu)*( P_0 - 1).*(P_0) + ((mu-1)/mu) - wbar + (1-rho)*beta*kappa.*Eta_1);

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

dY0_dA1 = (rho./kappa).*( (psi./mu).*( ( (lambda_1/delta).*(R_0.^(2/delta) )./Rbar).*( 2.*(R_0./Rbar) - (R_0./Rbar).^((delta-1)/delta) ) ) + (1-rho)*beta*kappa*theta_1 );

PS_Fraction = (rho./kappa).*( (psi./mu).*( ( (lambda_1/delta).*(R_0.^(2/delta) )./Rbar).*( 2.*(R_0./Rbar) - (R_0./Rbar).^((delta-1)/delta) ) ))./dY0_dA1;


%% Repeat when shutting off PS channel

chi=1;

theta_0 = -wbar/kappa;
theta_1 = (1-(1/mu))/kappa;

lambda_0 = beta*(1 + rho*(chi^(-gamma) - 1)*(1-theta_0));
lambda_1 = beta*rho*(chi^(-gamma) - 1)*theta_1;

Eta_1_Star = theta_0 + theta_1.*A_1;

M_1_Star = lambda_0 - lambda_1.*A_1;

R_0_Star = 1./M_1_Star;

Eta_0_Star = (1./kappa).*( (psi./mu)*( (R_0_Star./Rbar).^delta - 1).*(R_0_Star./Rbar).^delta + ((mu-1)/mu) - wbar + (1-rho)*beta*kappa.*Eta_1_Star);

N_0_Star = 1 - rho + rho.*Eta_0_Star;

Y_0_Star = N_0_Star;

V_0_Star = rho.*Eta_0_Star.^2;

P_0_Star = (R_0_Star/Rbar).^(1/delta);





%% Plot of Y_0 vs A_1

figure;
plot(100.*log(A_1), 100.*log(Y_0./Y_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(Y_0_Star./Y_0_Star(K)), 'color','blue','LineWidth',2,'LineStyle','-.');   % dashed blue line
xlabel('$\hat{A}_1$','interpreter','LaTeX','FontSize',15);
ylabel('$\hat{Y}_0$','interpreter','LaTeX','FontSize',15,'Rotation',0);
xlim([-10 10])
title('Output', 'interpreter', 'latex', 'FontSize', 16);
grid on

print('Figure1','-dpng');


%% Plots of other variables vs A_1


figure;

subplot(2,2,1)
plot(100.*log(A_1), 100.*log(Eta_0./Eta_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(Eta_0_Star./Eta_0_Star(K)),'color','blue','LineWidth',2,'LineStyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{\eta}_0$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
title('Job Finding Rate', 'interpreter', 'latex', 'FontSize', 16);
grid on

subplot(2,2,2)
plot(100.*log(A_1), 100.*log(V_0./V_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(V_0_Star./V_0_Star(K)),'color','blue','LineWidth',2,'LineStyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{v}_0$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
title('Vacancies', 'interpreter', 'latex', 'FontSize', 16);
grid on


subplot(2,2,3)
plot(100.*log(A_1), 100.*log(R_0./R_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(R_0_Star./R_0_Star(K)),'color','blue','LineWidth',2,'LineStyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{R}_{0}^{f}$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
title('Nominal Interest Rate', 'interpreter', 'latex', 'FontSize', 16);
grid on


subplot(2,2,4)
plot(100.*log(A_1), 100.*log(P_0./P_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(P_0_Star./P_0_Star(K)),'color','blue','LineWidth',2,'LineStyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{\Pi}_{0}$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
title('Inflation', 'interpreter', 'latex', 'FontSize', 16);
grid on

print('Figure2','-dpng');


%% Derivative calculations

A_1 = A_1';
Eta_0 = Eta_0';
Eta_0_Star = Eta_0_Star';
V_0 = V_0'
V_0_Star = V_0_Star';

b = regress(Eta_0, [ones(41,1), A_1 ] );

d_Eta_0 = b(2)

b = regress(Eta_0_Star, [ones(41,1), A_1 ] );

d_Eta_0_Star = b(2)

1 - d_Eta_0_Star/d_Eta_0

b = regress(V_0, [ones(41,1), A_1 ] );

d_V_0 = b(2)

b = regress(V_0_Star, [ones(41,1), A_1 ] );

d_V_0_Star = b(2)

1 - d_V_0_Star/d_V_0

%% Calculate MPC of employed agent

eta_bar = Mean_Eta1;

chi = 0.79;

rho_eps_bar = (1 - rho*(1-eta_bar)) + rho*(1-eta_bar)*(chi^-gamma);

R = 1/(beta*rho_eps_bar);

% Define the function to be solved
equation = @(x) x/(1-x) - (R/rho_eps_bar)*(rho*(1-eta_bar)*(chi^(-gamma-1)) + (1-rho*(1-eta_bar))*x)  ;

% Initial guess
initial_guess = 0.5;

% Solve the equation using fsolve
mpc = fsolve(equation, initial_guess);

disp(mpc);

average_mpc = Mean_Y0*mpc + (1-Mean_Y0)
