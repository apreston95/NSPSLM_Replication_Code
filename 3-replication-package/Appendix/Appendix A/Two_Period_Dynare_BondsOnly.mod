load params

// variables and parameters
var C_0_E Eta_0 Eta_1 C_1_EE C_1_EU  C_0_U C_1_UU C_1_UE q_0 R_0 P_0 P_1 Pi_1 B_0_E B_0_U N_0 N_1 e_0 e_1 Y_0 Y_1 MC_0 MC_1  junk;
parameters gamma rho beta wbar chi mu psi kappa A_1 Rbar;

// parameter values

gamma = 3;
rho = 0.05;
beta = 0.965;
wbar = 0.52;
mu = 6;
psi = 200;
kappa = 1;
Rbar = 1;
set_param_value('chi',par.chi);
set_param_value('A_1',par.A_1);

// Model declaration
model;

// Bond Euler employed

(C_0_E)^(-gamma) = beta*( (1-rho*(1-Eta_1(+1)))*(C_1_EE(+1))^(-gamma) + rho*(1-Eta_1(+1))*(C_1_EU(+1))^(-gamma) )*(R_0/Pi_1) ;

// Bond Euler unemployed

(C_0_U)^(-gamma) = beta*( (1-Eta_1(+1))*(C_1_UU(+1))^(-gamma) + Eta_1(+1)*(C_1_UE(+1))^(-gamma) )*(R_0/Pi_1) ;

// Employed Budget Constraint in 0

P_0*C_0_E = P_0*wbar - q_0*B_0_E;

// Unemployed Budget Constraint in 0

P_0*C_0_U = P_0*chi*wbar - q_0*B_0_U ;

// Employed Budget Constraint in 1 for employed in 0

P_1*C_1_EE = P_1*wbar + B_0_E ;

// Employed Budget Constraint in 1 for unemployed in 0

P_1*C_1_UE = P_1*wbar + B_0_U;

// Unmployed Budget Constraint in 1 for unemployed in 0

P_1*C_1_UU = P_1*chi*wbar + B_0_U;

// Unmployed Budget Constraint in 1 for employed in 0

P_1*C_1_EU = P_1*chi*wbar + B_0_E ;

// Asset Market Clearing

N_0*B_0_E + (1-N_0)*B_0_U = 0;

// Optimality of firms

mu*MC_0 = psi*(P_0-1)*P_0 + mu - 1;

MC_0 = wbar + kappa*Eta_0 - (1-rho)*beta*kappa*Eta_1;

mu*MC_1 = psi*(Pi_1-1)*Pi_1 + mu - 1;

MC_1 = (1/A_1)*(wbar + kappa*Eta_1);

// Employment LOM

N_0 = 1 - rho + Eta_0*e_0;

N_1 = (1-rho)*N_0 + Eta_1*e_1;

// Searcher LOM

e_0 = rho;

e_1 = 1 - N_0 + rho*N_0;

// Output

Y_0 = N_0;

Y_1 = A_1*N_1;


// Monetary Policy

R_0/Rbar = P_0;

Pi_1 = 1;

// Identities

q_0 = (1/R_0);

Pi_1 = P_1/P_0;


junk=0.9*junk(-1);

end;


initval;

N_0 = 0.95;
Y_0 = N_0;
e_0 = rho;
e_1 = 1 - N_0 + rho*N_0;
Eta_0 = 0.6;
Eta_1 = 0.6;
N_1 = (1-rho)*N_0 + Eta_1*e_1;
Y_1 = A_1*N_1;
MC_0 = 1;
MC_1 = 1;
B_0_E = 0.1;
B_0_U = -0.1;
P_0 = 1;
P_1 = 1;
Pi_1=1;
R_0 = 1;
C_1_EU = chi*wbar + B_0_E;
C_1_UU = chi*wbar + B_0_U;
C_1_UE = wbar + B_0_U ;
C_1_EE = wbar + B_0_E ;
q_0 = 0.99;
C_0_U = chi*wbar - q_0*B_0_U;
C_0_E = wbar - q_0*B_0_E;
junk=0;

end;

resid;

check;
steady(solve_algo=1, maxit=100000000);
// stoch_simul(order=3,irf=40, periods = 30000, pruning); 
 //stoch_simul(order=1,irf=40, periods = 1000); 