// News Shock Model W/ heterogeneity, written by Andy Preston  //

load params

@#define KBAR = 40



Var W R Pi Eta b Y A N M Eta_V MC U UN V JC C C_UU C_E_0 B_E_0 Psi_E_0 Ave_Cons_Loss logY logA logEta SavingRate URisk

@#for k in 1:KBAR

 C_E_@{k}  B_E_@{k} C_EU_@{k}  Psi_E_@{k}   

@#endfor

;


Varexo Eps_A Eps_A4 Eps_MPC;

Parameters 	Gamma
		Beta
		Rho
		Kappa
		Sigma
		Theta
		Delta_Pi
		Alpha
		Lambda
		Rho_A

        Rho_W
		
		W_Bar
		Pi_Bar
		Eta_Bar
		b_Bar
		Y_Bar
		A_Bar
		N_Bar
		M_Bar
		Eta_V_Bar
		MC_Bar
		U_Bar
		V_Bar
		B
		JC_Bar
		;

//****************************************************************************
//Set parameter values
//****************************************************************************

set_param_value('Gamma',par.Gamma);
set_param_value('Beta',par.Beta);
set_param_value('Rho',par.Rho);
set_param_value('Kappa',par.Kappa);
set_param_value('Sigma',par.Sigma);
set_param_value('Theta',par.Theta);
set_param_value('Delta_Pi',par.Delta_Pi);
set_param_value('Alpha',par.Alpha);
set_param_value('Lambda',par.Lambda);
set_param_value('Rho_A',par.Rho_A);

set_param_value('W_Bar',par.W_Bar);
set_param_value('Pi_Bar',par.Pi_Bar);
set_param_value('Eta_Bar',par.Eta_Bar);
set_param_value('b_Bar',par.b_Bar);
set_param_value('Y_Bar',par.Y_Bar);
set_param_value('A_Bar',par.A_Bar);
set_param_value('N_Bar',par.N_Bar);
set_param_value('M_Bar',par.M_Bar);
set_param_value('Eta_V_Bar',par.Eta_V_Bar);
set_param_value('MC_Bar',par.MC_Bar);
set_param_value('U_Bar',par.U_Bar);
set_param_value('V_Bar',par.V_Bar);
set_param_value('B',par.B);
set_param_value('JC_Bar',par.JC_Bar);

Rho_W = 0.9;


//****************************************************************************
//Model
//****************************************************************************

model;


//****************************************************************************
// HETEROGENEITY BLOCK
//****************************************************************************

// Employed BCs
C_E_0 + B_E_0 = W + Eps_MPC ;

@#for k in 1:KBAR-1

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*(B_E_@{k-1}(-1) + Eps_MPC) ;

@#endfor

@#for k in KBAR:KBAR

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*(B_E_@{k}(-1) + Eps_MPC) ;

@#endfor

// Employed Euler 
@#for k in 0:KBAR-1

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-Rho*(1-Eta(+1)))*(C_E_@{k+1}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k+1}(+1))^(-Gamma));

@#endfor

@#for k in KBAR:KBAR

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-Rho*(1-Eta(+1)))*(C_E_@{k}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k}(+1))^(-Gamma));

@#endfor

// Newly unemployed BC
@#for k in 1:KBAR

C_EU_@{k} = (R(-1)/(1+Pi))*(B_E_@{k-1}(-1)+ Eps_MPC) + b  ;

@#endfor

// Continuing unemployed BC

C_UU = b;




// Population shares

Psi_E_0 = M;

@#for k in 1:KBAR-1
 log(Psi_E_@{k}) = log(M(-@{k})) 
                @#for i in 1:k 
                 + log(JC(-@{i}))
                 @#endfor
                 ;

@#endfor




Psi_E_@{KBAR} = N  
                    @#for i in 0:KBAR-1 
                  - Psi_E_@{i}
                  @#endfor
                  ;



 
// Bond market clearing

B = 1*(
       @#for k in 0:KBAR 
        + Psi_E_@{k}*B_E_@{k}
          @#endfor
                  );


// Average cons. loss

Ave_Cons_Loss = (Psi_E_@{KBAR}/N)*(C_EU_@{KBAR}/C_E_@{KBAR}) +

            1*(
            @#for k in 1:KBAR 
        + (Psi_E_@{k-1}/N)*(C_EU_@{k}/C_E_@{k-1})
          @#endfor
                    );

// Average saving rate

SavingRate = (Psi_E_@{KBAR}/N)*(1-((C_E_@{KBAR})/(W+(R(-1)/(1+Pi))*B_E_@{KBAR}(-1)))) + 

            1*(
            @#for k in 1:KBAR-1 
        + (Psi_E_@{k}/N)*(1-((C_E_@{k})/(W+(R(-1)/(1+Pi))*B_E_@{k-1}(-1))))
          @#endfor
                    )
            + (Psi_E_0/N)*(1-(C_E_0/W))
                        
                       ;
		 
                  
//****************************************************************************
// NON-HETEROGENEITY BLOCK
//****************************************************************************

URisk = 1 - ((1-Rho*(1-Eta(+1)))*(1-Rho*(1-Eta(+2)))*(1-Rho*(1-Eta(+3)))*(1-Rho*(1-Eta(+4))));

JC = 1 - Rho*(1-Eta(+1));

[name='Prod. Function']
Y = A*N;

[name='Employment LOM']
N = (1-Rho)*N(-1) + M;

[name='Unmployment LOM']
U = 1 - (1-Rho)*N(-1);

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
R/steady_state(R) =((1+Pi)/(1+Pi_Bar))^(Delta_Pi);

[name='Matching Function']
M = U^(Alpha)*V^(1-Alpha);

[name='Agg. Consumption']
C = Y - Kappa*V;

[name='Home Production']
b = Lambda*steady_state(W);

[name='TFP']
log(A) = Rho_A*log(A(-1)) + Eps_A + Eps_A4(-4);

logY = log(Y);
logA = log(A);
logEta = log(Eta);
UN = 1-N;


end; 




initval(all_values_required);

W = W_Bar;
Pi = Pi_Bar;
Eta = Eta_Bar;
b = b_Bar;
Y = Y_Bar;
A = A_Bar;
logY = log(Y);
logA = log(A);
N = N_Bar;
M = M_Bar;
Eta_V = Eta_V_Bar;
MC = MC_Bar;
U = U_Bar;
UN = 1-N;
V = V_Bar;
C = Y - Kappa*V_Bar;
JC = JC_Bar;
Eps_A = 0;
Eps_A4 = 0;
Eps_MPC =0;

Psi_E_0 = M;

@#for k in 1:KBAR-1
 Psi_E_@{k} = M*(JC)^(@{k});
                 
@#endfor


@#for k in KBAR:KBAR

Psi_E_@{k} = N
                  @#for i in 0:KBAR-1 
                  - Psi_E_@{i}
                  @#endfor
                  ;
@#endfor
		 
		 
R = 1/Beta;

C_UU = b;

@#for k in 0:KBAR
C_E_@{k} = W;
@#endfor

@#for k in 1:KBAR
C_EU_@{k} = 0.8*C_E_@{k};
@#endfor

B_E_0 = 0;

@#for k in 1:KBAR
B_E_@{k} = B;
@#endfor

Ave_Cons_Loss = 0.8;

SavingRate = 0.1;

URisk = 1 - ((1-Rho*(1-Eta))*(1-Rho*(1-Eta))*(1-Rho*(1-Eta))*(1-Rho*(1-Eta)));

logEta = log(Eta);

end;


shocks;
var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;
var Eps_MPC; stderr 0;
end;

resid;

steady(solve_algo=4, maxit=100000000);

model_diagnostics;



stoch_simul(order=3, periods=0, irf=0); % no default IRFs, we'll compute them manually






