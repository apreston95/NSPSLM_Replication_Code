clc
clear all
close all

KBAR = 80;

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);


Gamma = 2;
Beta = 0.99;
Rho = 0.044;
Sigma = 6;
Theta = 58.7;
Delta_Pi = 1.5;
Alpha = 0.65;
Lambda = 0.717;
Rho_A = 0.99;
Vacfrac = 0.01;
Delta_News = 0;



A_Bar = 1;
Pi_Bar = 0;
N_Bar = 0.955;
U_Bar = 1-N_Bar*(1-Rho);
M_Bar = Rho*N_Bar;
Eta_Bar = M_Bar/U_Bar;
V_Bar = (M_Bar/U_Bar^(Alpha))^(1/(1-Alpha));
Eta_V_Bar = M_Bar/V_Bar;
Y_Bar = A_Bar*N_Bar;
B = 0.065*Y_Bar;
MC_Bar = (Sigma-1)/Sigma;
Kappa = Vacfrac*(Y_Bar/V_Bar);
Wage = @(x) MC_Bar - (1/A_Bar)*(x + (Kappa/Eta_V_Bar) - Beta*(1-Rho)*(Kappa/Eta_V_Bar));
W_Bar = fsolve(Wage,1,options);
b_Bar = Lambda*W_Bar;
JC_Bar = 1 - Rho*(1-Eta_Bar);

Psi_E_vec = zeros(1,KBAR);
for j = 0:KBAR-2
  Psi_E_vec(j+1) = M_Bar*(JC_Bar)^(j);
end


Psi_E_vec(1,KBAR) = N_Bar - sum(Psi_E_vec(1,1:KBAR-1));


par.Gamma = Gamma;
par.Beta = Beta;
par.Rho = Rho;
par.Sigma = Sigma;
par.Theta = Theta;
par.Delta_Pi = Delta_Pi;
par.Alpha = Alpha;
par.Lambda = Lambda;
par.Rho_A = Rho_A;
par.Kappa = Kappa;

par.A_Bar = A_Bar;
par.Pi_Bar = Pi_Bar;
par.N_Bar = N_Bar;
par.U_Bar = U_Bar;
par.M_Bar = M_Bar;
par.Eta_Bar = Eta_Bar;
par.Y_Bar = Y_Bar;
par.B = B;
par.MC_Bar = MC_Bar;
par.V_Bar = V_Bar;
par.Eta_V_Bar = Eta_V_Bar;
par.W_Bar = W_Bar;
par.b_Bar = b_Bar;
par.JC_Bar = JC_Bar;
par.Delta_News = Delta_News;



save params par;

%% Baseline (Figure 4)


dynare NSPSLM_Het noclearall nolog;



%% Countercyclical Unemployment Benefit and More Hawksish Monetary Policy (Figure 6)

dynare NSPSLM_Het_OptimalPolicy noclearall  nolog;



