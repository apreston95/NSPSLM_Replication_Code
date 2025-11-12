clc
clear all
close all

KBAR = 80;

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);


Gamma = 2;
Beta = 0.994;
PEU = 0.0228;
Sigma = 6;
Theta = 58.7;
Delta_Pi = 1.5;
Lambda = 0.735;
Rho_A = 0.99;
Kappa1 = 1;


A_Bar = 1;
Pi_Bar = 0;
U_Bar = 0.045;
PUE = PEU*(1/U_Bar - 1);
MC_Bar = (Sigma - 1)/Sigma;
W_Bar = A_Bar*MC_Bar;
b_Bar = Lambda*W_Bar*0.33;

Psi_E_Vec = zeros(1,KBAR+1);

for j = 0:KBAR-1
Psi_E_Vec(j+1) = PUE*((1-PEU)^(j+1))*U_Bar
end

Psi_E_Vec(KBAR+1) = 1 - U_Bar - sum(Psi_E_Vec(1:KBAR));
%% 



par.Gamma = Gamma;
par.Beta = Beta;
par.PEU = PEU;
par.Sigma = Sigma;
par.Theta = Theta;
par.Delta_Pi = Delta_Pi;
par.Lambda = Lambda;
par.Rho_A = Rho_A;
par.PUE = PUE;
par.Kappa1 = Kappa1;


par.A_Bar = A_Bar;
par.Pi_Bar = Pi_Bar;
par.U_Bar = U_Bar;
par.MC_Bar = MC_Bar;
par.W_Bar = W_Bar;
par.b_Bar = b_Bar;

save params par;

%% 


dynare NSPSLM_Spot_Het noclearall;


