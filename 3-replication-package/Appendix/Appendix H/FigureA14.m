clear all
close all

T=100000;

kappa = 1;
mu = 6;
rho = 0.05;
chi = 0.79;
gamma = 3;
wbar = 0.52;
beta = 0.965;
beta_A = 1;


sigma=0.5;

theta_0 = -wbar/kappa;
theta_1 = (1-(1/mu))/kappa;


tol = 0.00001;

A_Grid = 0.9:0.005:1.1;

K = (length(A_Grid)+1)/2 ;

for i = 1:length(A_Grid)
    
A_1_tilde = A_Grid(i);

Eta_1 = theta_0 + theta_1*A_1_tilde;

eta_bar = Eta_1;

Eta_0 = (1/kappa)*( (mu-1)/mu - wbar + (1-rho)*beta*kappa*Eta_1 );

N_0 = 1 - rho + rho*Eta_0;

w_e = 0.0001;

w_u = -w_e*(N_0/(1-N_0));

diff = 1;

R_0 = 1;

Q_0 = 0.8;

Cov = (1/Q_0)*(beta_A*theta_1*sigma^2);

E_R_SQ = (1/Q_0^2)*(1 + (beta_A^2)*(sigma^2) ); 

iter = 0;

while diff > tol


c_ee = wbar + R_0*w_e;
c_eu = chi*wbar + R_0*w_e;

c_uu = chi*wbar + R_0*w_u;
c_ue = wbar + R_0*w_u;

chi_tilde = c_eu/c_ee;
chi_hat = c_uu/c_ue;

delta_e = c_ee/( gamma * E_R_SQ * (1 + rho*(1-eta_bar)*( chi_tilde^(-gamma-1) - 1) ));

phi_e_r_eta = rho*(1 - chi_tilde^(-gamma) );

phi_e_r = 1 + rho*(1-eta_bar)*(chi_tilde^(-gamma) - 1);

delta_u = c_eu/( gamma * E_R_SQ * (eta_bar + (1-eta_bar)*(chi_hat)^(-gamma-1) ) );

phi_u_r_eta = 1 - chi_hat^(-gamma) ;

phi_u_r = eta_bar + (1-eta_bar)*(chi_hat)^(-gamma) ;

EP = -Cov* (  (N_0*delta_e*phi_e_r_eta + (1-N_0)*delta_u*phi_u_r_eta)/( N_0*delta_e*phi_e_r + (1-N_0)*delta_u*phi_u_r ) )

%% Portfolios


a_0_e_r = delta_e*(phi_e_r_eta*Cov + phi_e_r*EP)

alpha_0_e = a_0_e_r/w_e

a_0_e_f = (1-alpha_0_e)*w_e

a_0_u_r = delta_u*(phi_u_r_eta*Cov + phi_u_r*EP)

alpha_0_u = a_0_u_r/w_u

a_0_u_f = (1-alpha_0_u)*w_u

excess_r = N_0*a_0_e_r + (1-N_0)*a_0_u_r

excess_f = N_0*a_0_e_f + (1-N_0)*a_0_u_f

MPR_e = (1/( gamma * E_R_SQ * (1 + rho*(1-eta_bar)*( chi_tilde^(-gamma-1) - 1) )))*(phi_e_r_eta*Cov + phi_e_r*EP)

MPR_u = (1/( gamma * E_R_SQ * (eta_bar + (1-eta_bar)*(chi_hat)^(-gamma-1) ) ))*(phi_u_r_eta*Cov + phi_u_r*EP)

c_0_e = wbar - w_e;


c_0_u = chi*wbar - w_u;

%% Solve for prices and moments to verify guesses


q_0 = c_0_e^(gamma)*beta*( (1 - rho*(1-eta_bar))*(c_ee)^(-gamma) + rho*(1-eta_bar)*(c_eu)^(-gamma) - gamma*(1 - rho*(1-eta_bar))*(c_ee)^(-gamma-1)*w_e*alpha_0_e*EP - gamma*rho*(1-eta_bar)*(c_eu)^(-gamma-1)*w_e*alpha_0_e*EP );

R_0_sol = 1/q_0;

E_R_1 = EP + R_0_sol;

Q_0 = 1/E_R_1;

Cov_sol = (1/Q_0)*(beta_A*theta_1*sigma^2);

E_R_SQ_sol = (1/Q_0^2)*(1 + (beta_A^2)*(sigma^2) );

diff = abs(R_0_sol - R_0)

Cov = Cov_sol;

R_0 = R_0_sol;

E_R_SQ = E_R_SQ_sol;

%w_e = a_0_e_r + a_0_e_f;

equation = @(w_u_sol) q_0 - (chi*wbar - w_u_sol)^gamma * beta * ((1-eta_bar)*(chi*wbar + R_0*w_u_sol)^(-gamma) + eta_bar*(wbar + R_0*w_u_sol)^(-gamma) - gamma*(1-eta_bar)*(chi*wbar + R_0*w_u_sol)^(-gamma-1)*w_u_sol*alpha_0_u*EP - gamma*eta_bar*(wbar + R_0*w_u_sol)^(-gamma-1)*w_u_sol*alpha_0_u*EP);

% Initial guess for fsolve
initial_guess = w_u;

% Solve the equation using fsolve
w_u = fsolve(equation, initial_guess);

w_e = -((1-N_0)/N_0)*w_u;

%diff = abs(excess_w);

iter = iter + 1;

end;



EP_Grid(i) = EP;

PVS_Grid(i) = log(Q_0/q_0);

a_e_r_Grid(i) = a_0_e_r;

a_u_r_Grid(i) = a_0_u_r;

alpha_e_grid(i) = a_0_e_r/(a_0_e_r + a_0_e_f);

alpha_u_grid(i) = a_0_u_r/(a_0_u_r + a_0_u_f);

R_0_Grid(i) = R_0;

excess_r_grid(i) = excess_r;

excess_f_grid(i) = excess_f;

iter_grid(i) = iter;

w_e_grid(i) = w_e;

w_u_grid(i) = w_u;

eta_grid(i) = eta_bar;

N_0_grid(i) = N_0;



end;

%% Figures


figure;
plot(A_Grid,400*EP_Grid,'color','red','LineWidth',2);
hold on
xlim([0.9 1.1])
xlabel('$\widetilde{A}_{1}$','interpreter','LaTeX', 'FontSize', 15);
ylabel('Annualised Equity Premium (\%)','interpreter','LaTeX', 'FontSize', 15);
% 
% print('ER_PC','-dpng');

PVS_SS = PVS_Grid(K);

figure;
% subplot(1,2,1)
plot(A_Grid,100*(PVS_Grid ),'color','red','LineWidth',2.5);
hold on
xlim([0.9 1.1])
xlabel('$\widetilde{A}_{1}$','interpreter','LaTeX', 'FontSize', 15);
ylabel('PVS (\%)','interpreter','LaTeX', 'FontSize', 15);
title('PVS','interpreter','LaTeX', 'FontSize', 18);

% subplot(1,2,2)
% plot(A_Grid,a_e_r_Grid,'color','blue','LineWidth',2.5);
% hold on
% plot(A_Grid,a_u_r_Grid,'color','green','LineWidth',2.5);
% hold on
% xlim([0.9 1.1])
% xlabel('$\widetilde{A}_{1}$','interpreter','LaTeX', 'FontSize', 15);
% ylabel('$a_0^{i,r}$','interpreter','LaTeX', 'FontSize', 15);
% legend('$a_0^{E,r}$','$a_0^{U,r}$')


print('FigureA14','-dpng');

% figure;
% plot(A_Grid,1200*(R_0_Grid-1),'color','blue','LineWidth',2);
% hold on
% xlim([0.9 1.1])
% xlabel('$\widetilde{A}_{1}$ ','interpreter','LaTeX', 'FontSize', 15);
% ylabel('Annualised Risk-Free Rate (\%)','interpreter','LaTeX', 'FontSize', 15);
% 
% print('RF_PC','-dpng');
