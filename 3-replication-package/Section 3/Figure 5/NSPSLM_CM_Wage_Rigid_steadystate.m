function [ys,params,check] = NSPSLM_2021_CM_Wage_Rigid(ys,exo,M_,options_)
% function [ys,params,check] = NSPSLM_2021(ys,exo,M_,options_)
% computes the steady state for the NSPSLM_2021.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);

A = A_Bar;
Pi = Pi_Bar;
N = N_Bar;
U = 1 - (1-Rho)*N;
M = Rho*N;
Eta = M/U;
MC = (Sigma-1)/Sigma;
Y = A*N;
V = (M/U^(Alpha))^(1/(1-Alpha));
Kappa = Vacfrac*(Y/V);
Eta_V = M/V;
Wage = @(x) MC - (1/A)*(x + (Kappa/Eta_V) - Beta*(1-Rho)*(Kappa/Eta_V));
W = fsolve(Wage,1,options);

b = Lambda*W;
C  = Y - ( Kappa*V + (Theta/2)*(Pi^2)*Y);


R = 1/Beta;

logC = log(C);
logA = log(A);
logY = log(Y);
Real = R/(1+Pi);
Urisk = 1 - ((1-Rho*(1-Eta))*(1-Rho*(1-Eta))*(1-Rho*(1-Eta))*(1-Rho*(1-Eta)));


R_Bar = R;
W_Bar = W;
Y_Bar = Y;

%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end





