
load params

@#define KBAR = 80

Var W R Pi B b Y A N MC U logY logA Ave_Cons_Loss C_UU C_E_0 B_E_0 N_0 Psi_E_0 
        
@#for k in 1:KBAR

 C_E_@{k}  B_E_@{k} C_EU_@{k}  N_@{k} Psi_E_@{k}  

@#endfor

;


Varexo Eps_A Eps_A4;

Parameters 	Gamma
		Beta
		PEU
        PUE
		Kappa0
        Kappa1
		Sigma
		Theta
		Delta_Pi
		Lambda
		Rho_A
        Bfrac

        U_Bar
		W_Bar
		Pi_Bar
   		A_Bar
        MC_Bar
        b_Bar


		;

//****************************************************************************
//Set parameter values
//****************************************************************************

set_param_value('Gamma',par.Gamma);
set_param_value('Beta',par.Beta);
set_param_value('PEU',par.PEU);
set_param_value('PUE',par.PUE);
set_param_value('Sigma',par.Sigma);
set_param_value('Theta',par.Theta);
set_param_value('Delta_Pi',par.Delta_Pi);
set_param_value('Lambda',par.Lambda);
set_param_value('Rho_A',par.Rho_A);
set_param_value('Kappa1',par.Kappa1);


set_param_value('W_Bar',par.W_Bar);
set_param_value('Pi_Bar',par.Pi_Bar);
set_param_value('b_Bar',par.b_Bar);
set_param_value('A_Bar',par.A_Bar);
set_param_value('MC_Bar',par.MC_Bar);
set_param_value('U_Bar',par.U_Bar);

Kappa0 = 27.8;
Bfrac = 0.08;

//****************************************************************************
//Model
//****************************************************************************

model;

//****************************************************************************
//Heterogeneity Block
//****************************************************************************

// Employed BCs
C_E_0 + B_E_0 = W*N_0 ;

@#for k in 1:KBAR-1

C_E_@{k} + B_E_@{k} = W*N_@{k} + R(-1)/(1+Pi)*B_E_@{k-1}(-1) ;

@#endfor

@#for k in KBAR:KBAR

C_E_@{k} + B_E_@{k} = W*N_@{k} + R(-1)/(1+Pi)*B_E_@{k}(-1);

@#endfor

// Employed Euler 
@#for k in 0:KBAR-1

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-PEU)*(C_E_@{k+1}(+1))^(-Gamma) + PEU*(C_EU_@{k+1}(+1))^(-Gamma));

@#endfor

@#for k in KBAR:KBAR

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-PEU)*(C_E_@{k}(+1))^(-Gamma) + PEU*(C_EU_@{k}(+1))^(-Gamma));

@#endfor

// Employed Labour Supply

@#for k in 0:KBAR

W*(C_E_@{k})^(-Gamma) = Kappa0*(N_@{k})^(Kappa1);

@#endfor

// Newly unemployed BC
@#for k in 1:KBAR

C_EU_@{k} = (R(-1)/(1+Pi))*B_E_@{k-1}(-1) + b;

@#endfor

// Continuing unemployed BC

C_UU = b;


// Population Shares


@#for k in 0:KBAR-1

Psi_E_@{k} = PUE*(1-PEU)^(@{k})*U;

@#endfor

Psi_E_@{KBAR} = 1 - U - 1*(
                        @#for k in 0:KBAR-1 
                        + Psi_E_@{k}
                        @#endfor
                                             );


// Bond market clearing

B = Bfrac*steady_state(Y);

B = 1*(
       @#for k in 0:KBAR 
        + Psi_E_@{k}*B_E_@{k}
          @#endfor
                  );

// Labour market clearing

N = 1*(
       @#for k in 0:KBAR 
        + Psi_E_@{k}*N_@{k}
          @#endfor
                  );


// Average cons. loss

Ave_Cons_Loss = (Psi_E_@{KBAR}/(1-U))*(C_EU_@{KBAR}/C_E_@{KBAR}) +

            1*(
            @#for k in 1:KBAR 
        + (Psi_E_@{k-1}/(1-U))*(C_EU_@{k}/C_E_@{k-1})
          @#endfor
                    );

         

//****************************************************************************
//Non-Heterogeneity Block
//****************************************************************************


[name='MC']
MC = W/A;

[name='Phillips Curve']
1 - Sigma + Sigma*MC = Theta*(Pi+1)*Pi - Theta*Beta*((Pi(+1)+1)*Pi(+1)*(Y(+1)/Y));

[name='Prod. Function']
Y = A*N;

[name='Unemployment']
U = PEU/(PEU+PUE);

[name='Monetary Policy']
R/steady_state(R) =((1+Pi)/(1+Pi_Bar))^(Delta_Pi);

[name='Home Production']
b = Lambda*steady_state(W)*steady_state(N);

[name='TFP']
log(A) = Rho_A*log(A(-1)) - Eps_A - Eps_A4(-4);

logA = log(A);
logY = log(Y);


end; 




initval(all_values_required);

W = W_Bar;
Pi = Pi_Bar;
b = b_Bar;
Y = 0.33;
B = Bfrac*Y;
N = 0.33;
A = A_Bar;
logY = log(Y);
logA = log(A);
MC = MC_Bar;
U = U_Bar;
Eps_A = 0;
Eps_A4 = 0;


@#for k in 0:KBAR-1
 Psi_E_@{k} = PUE*(1-PEU)^(@{k})*U;
                 
@#endfor


Psi_E_@{KBAR} = 1 - U - 1*(
                        @#for k in 0:KBAR-1 
                        + Psi_E_@{k}
                        @#endfor
                                             );
		 
		 
R = 1/Beta;

C_UU = b;

@#for k in 0:KBAR
C_E_@{k} = 0.7*Y;
@#endfor

@#for k in 0:KBAR
N_@{k} = N;
@#endfor

@#for k in 1:KBAR
C_EU_@{k} = Lambda*C_E_@{k};
@#endfor

@#for k in 0:KBAR
B_E_@{k} = B;
@#endfor

Ave_Cons_Loss = 0.8;

end;


shocks;
var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;
end;

resid;

steady(solve_algo=1, maxit=100000000);

model_diagnostics;

stoch_simul(order=1, nocorr,nograph, nomoments,irf=21);

Hor = 20;

IRF_TFP_NEWS4 = oo_.irfs.logA_Eps_A4;
IRF_OUTPUT_NEWS4 = oo_.irfs.logY_Eps_A4;
IRF_NOMINAL_NEWS4 = oo_.irfs.R_Eps_A4;
IRF_INFLATION_NEWS4 = oo_.irfs.Pi_Eps_A4;



figure;


subplot(2,2,1)
plot(0:Hor,IRF_TFP_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 Hor],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('TFP','interpreter','LaTeX');


subplot(2,2,2)
plot(0:Hor,IRF_OUTPUT_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 Hor],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Output','interpreter','LaTeX');



subplot(2,2,3)
plot(0:Hor,IRF_NOMINAL_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 Hor],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Nominal interest rate','interpreter','LaTeX');


subplot(2,2,4)
plot(0:Hor,IRF_INFLATION_NEWS4(1:end)*100,'color','black','LineWidth',2);
hold on
line([0 Hor],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Inflation','interpreter','LaTeX');

IRF_OUTPUT_NEWS4_SPOT_HET = oo_.irfs.logY_Eps_A4;
IRF_NOMINAL_NEWS4_SPOT_HET = oo_.irfs.R_Eps_A4;
save('IRFs_SPOT_HET.mat','IRF_OUTPUT_NEWS4_SPOT_HET','IRF_NOMINAL_NEWS4_SPOT_HET');
