
Var W R Pi Eta b Y A N M Eta_V MC U V C logC logY logA Real Urisk;

Varexo Eps_A Eps_A4;

Parameters 	Gamma
		Beta
		Rho
		Kappa
		Sigma
		Theta
		Chi
		Delta_Pi
		Alpha
		Lambda
		Rho_A
        Rho_W
        Phi
        Vacfrac

		N_Bar
		W_Bar
		R_Bar
		Pi_Bar
   		A_Bar
        Y_Bar
		;

//****************************************************************************
//Set parameter values
//****************************************************************************

Gamma = 2;
Beta = 0.99;
Rho = 0.044;
Sigma = 6;
Theta = 58.7;
Chi = 0.04;
Delta_Pi = 1.5;
Alpha = 0.65;
Lambda = 0.82;
Rho_A = 0.99;
Rho_W = 0.9;
Vacfrac = 0.01;

A_Bar = 1;
N_Bar = 0.955;
Pi_Bar = 0;

//****************************************************************************
//Model
//****************************************************************************
model;

[name='Euler equation']
C^(-Gamma) = Beta*(R/(1+Pi(+1)))*C(+1)^(-Gamma);

[name='Prod. Function']
Y = A*N;

[name='Employment LOM']
N = (1-Rho)*N(-1) + M;

[name='MC']
MC = (1/A)*(W + (Kappa/Eta_V) - Beta*(1-Rho)*(Kappa/Eta_V(+1)));

[name='Phillips Curve']
1 - Sigma + Sigma*MC = Theta*(Pi+1)*Pi - Theta*Beta*((Pi(+1)+1)*Pi(+1)*(Y(+1)/Y));

[name='Wage Setting']
log(W/W_Bar) = Rho_W*log(W(-1)/W_Bar) + (1-Rho_W)*log(A/A_Bar);

[name='Job Finding Rate']
Eta = M/U;

[name='Vacancy Filling Rate']
Eta_V = M/V;

[name='Monetary Policy']
log(R/R_Bar) = Delta_Pi*log((1+Pi)/(1+Pi_Bar));

[name='Matching Function']
M = U^(Alpha)*V^(1-Alpha);

[name='Unemployment']
U = 1 - (1-Rho)*N(-1);

[name='Home Production']
b = Lambda*steady_state(W);

[name='Goods Market Clearing']
C + Kappa*V + (Theta/2)*(Pi^2)*Y = Y;


[name='TFP']
log(A) = Rho_A*log(A(-1)) - Eps_A - Eps_A4(-4);

logA = log(A);
logY = log(Y);
logC = log(C);
Real = R/(1+Pi(+1));
Urisk = 1 - ((1-Rho*(1-Eta(+1)))*(1-Rho*(1-Eta(+2)))*(1-Rho*(1-Eta(+3)))*(1-Rho*(1-Eta(+4))));


end; 



shocks;

var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;

end;

steady;
check;

model_diagnostics;

stoch_simul(order=1, nocorr,nograph, nomoments,irf=21);

H = 20;

IRF_TFP_NEWS4 = oo_.irfs.logA_Eps_A4;
IRF_OUTPUT_NEWS4 = oo_.irfs.logY_Eps_A4;
IRF_CONSUMPTION_NEWS4 = oo_.irfs.logC_Eps_A4;
IRF_UNEMPLOYMENT_NEWS4 = oo_.irfs.U_Eps_A4;
IRF_VACANCIES_NEWS4 = oo_.irfs.V_Eps_A4;
IRF_NOMINAL_NEWS4 = oo_.irfs.R_Eps_A4;
IRF_REAL_NEWS4 = oo_.irfs.Real_Eps_A4;
IRF_INFLATION_NEWS4 = oo_.irfs.Pi_Eps_A4;
IRF_UR_NEWS4 = oo_.irfs.Urisk_Eps_A4;



figure;


subplot(2,4,1)
plot(0:H,IRF_TFP_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('TFP','interpreter','LaTeX');


subplot(2,4,2)
plot(0:H,IRF_OUTPUT_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Output','interpreter','LaTeX');

subplot(2,4,3)
plot(0:H,IRF_CONSUMPTION_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Consumption','interpreter','LaTeX');

subplot(2,4,4)
plot(0:H,IRF_UNEMPLOYMENT_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Unemployment','interpreter','LaTeX');

subplot(2,4,5)
plot(0:H,IRF_VACANCIES_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Vacancies','interpreter','LaTeX');

subplot(2,4,7)
plot(0:H,IRF_NOMINAL_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Nominal interest rate','interpreter','LaTeX');

subplot(2,4,6)
plot(0:H,IRF_UR_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Unemployment risk','interpreter','LaTeX');

subplot(2,4,8)
plot(0:H,IRF_INFLATION_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Inflation','interpreter','LaTeX');

IRF_OUTPUT_NEWS4_CM_RIGID = oo_.irfs.logY_Eps_A4;
IRF_NOMINAL_NEWS4_CM_RIGID = oo_.irfs.R_Eps_A4;
IRF_VACANCIES_NEWS4_CM_RIGID = oo_.irfs.V_Eps_A4;
save('IRFs_CM_RIGID.mat','IRF_VACANCIES_NEWS4_CM_RIGID','IRF_OUTPUT_NEWS4_CM_RIGID','IRF_NOMINAL_NEWS4_CM_RIGID');
