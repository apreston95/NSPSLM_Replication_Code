clear all
close all
clc

%% Parameters

kappa = 1;
mu = 6;
rho = 0.05;
chi = 0.79;
gamma = 3;
psi = 200;
wbar = 0.52;
beta = 0.965;
Rbar = 1;

alpha_0 = beta*(1 + rho*(chi^(-gamma) - 1 ) );
alpha_1 = beta*rho*(chi^(-gamma) - 1 );

theta_0 = -wbar*alpha_0/( kappa*(1 + wbar*alpha_1));
theta_1 = (1-1/mu)/( kappa*(1 + wbar*alpha_1));

lambda_0 = alpha_0 - alpha_1*theta_0;
lambda_1 = alpha_1*theta_1;


%% Variables


A_1 = 0.9:0.005:1.1;

Eta_1 = theta_0 + theta_1.*A_1;

M_1 = lambda_0 - lambda_1.*A_1;

R_0 = 1./M_1;

Eta_0 = (1./kappa).*( ((mu-1)/mu) - wbar*lambda_0 + (psi/mu).*(R_0-1).*R_0 + (1-rho)*beta*kappa*theta_0 + lambda_1*A_1.*wbar + (1-rho)*beta*kappa.*theta_1*A_1);

N_0 = 1 - rho + rho.*Eta_0;

Y_0 = N_0;


V_0 = rho.*Eta_0.^2;

V_1 = rho.*Eta_1.^2;

e_1 = 1 - N_0 + rho*N_0;

N_1 = (1-rho).* N_0 + Eta_1.*e_1;

K = 0.5*(length(A_1) - 1) + 1;

Mean_Eta1 = Eta_1(K);

Mean_Eta0 = Eta_0(K);

Mean_Y0 = Y_0(K);

Mean_R0 = R_0(K);

Mean_N1 = N_1(K);

dY0_dA1 = (rho./kappa).*( (psi/mu)*lambda_1.*(R_0).^2.*(2*R_0-1) +  lambda_1.*wbar + (1-rho)*beta*kappa.*theta_1 );

PS_Fraction = ((rho./kappa).*lambda_1)./dY0_dA1;

P_0 = R_0;








%% Repeat when shutting off PS channel

chi=1;

alpha_0 = beta*(1 + rho*(chi^(-gamma) - 1 ) );
alpha_1 = beta*rho*(chi^(-gamma) - 1 );

theta_0 = -wbar*alpha_0/( kappa*(1 + wbar*alpha_1));
theta_1 = (1-1/mu)/( kappa*(1 + wbar*alpha_1));

lambda_0 = alpha_0 - alpha_1*theta_0;
lambda_1 = alpha_1*theta_1;

Eta_1_Star = theta_0 + theta_1.*A_1;

M_1_Star = lambda_0 - lambda_1.*A_1;

R_0_Star = 1./M_1_Star;

Eta_0_Star = (1./kappa).*( ((mu-1)/mu) - wbar*lambda_0 + (psi/mu).*(R_0_Star-1).*R_0_Star + (1-rho)*beta*kappa*theta_0 + lambda_1*A_1.*wbar + (1-rho)*beta*kappa.*theta_1*A_1);

N_0_Star = 1 - rho + rho.*Eta_0_Star;

Y_0_Star = N_0_Star;

V_0_Star = rho.*Eta_0_Star.^2;

P_0_Star = R_0_Star/Rbar;



%% Plot of Y_0 vs A_1

figure;
%plot(log(A_1), log(Y_0/Mean_Y0),'color','black','LineWidth',2);
plot(100.*log(A_1), 100.*log(Y_0./Y_0(K)),'color','black','LineWidth',2);
hold on
plot(100.*log(A_1), 100*log(Y_0_Star./Y_0_Star(K)),'color','blue','LineWidth',2,'linestyle','-.');
xlabel('$\hat{A}_1$ ','interpreter','LaTeX', 'FontSize', 15);
ylabel('$\hat{Y}_0$','interpreter','LaTeX', 'FontSize', 15, 'Rotation', 0);
xlim( [-10 10])
grid on

print('FigureA2','-dpng');

