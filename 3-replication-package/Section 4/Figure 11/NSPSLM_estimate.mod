%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:
%   1) Compute the (nonlinear) steady state (user-provided via steady_state_model).
%   2) Use the log-linear model equations.
%   3) Estimate parameters (alpha, gamma, mu, psi, chi, kappa, delta, phi_A, rho_W)
%      plus the standard deviations of shocks eA, eA4.
%   4) Data: Observables are dy, r, pi, n.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% 1. Declare variables and shocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var w    // log-deviation of the real wage
    pi   // inflation
    r    // nominal interest rate
    y    // output
    n    // employment
    e    // number of job searchers
    v    // vacancies
    m    // matches
    f    // job filling rate
    eta  // job finding rate
    mc   // marginal cost
    a   // TFP
    u_R // Monetary policy disturbance
    u_Pi // Cost push disturbance
    dy  // Output growth as measurement equation
    u;   // Unemployment rate as potential measurement equation

varexo eA  ${\varepsilon^A}$ // TFP shock this period
       eA4 ${\varepsilon^{A,N_4}}$ // TFP shock from 4 periods ago
       eA8 ${\varepsilon^{A,N_8}}$  // TFP shock from 8 periods ago
       eR   ${\varepsilon^R}$ // Monetary shock 
       ePi ${\varepsilon^{\Pi}} $ ;   // Cost push shock


%%%%%%%%%%%% 2. Declare parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters alpha    ${\alpha}$ %: matching elasticity
           gamma    ${\gamma}$ % inverse IES
           beta     ${\beta}$ % discount factor (fixed below)
           Nbar     ${\overline{N}}$ % steady-state employment
           kappa    ${\kappa}$ % vacancy-posting cost
           mu       ${\mu}$ % markup-related parameter
           rho      ${\rho}$ % separation rate (fixed below)
           chi      ${\chi}$ % utility/wage parameter
           psi      ${\psi}$ % NKPC slope parameter
           rho_W    ${\rho_W}$ % wage rule inertia
           rho_A    ${\rho_A}$ % wage rule loading on TFP
           delta    ${\delta}$ % Taylor rule coefficient
           phi_A    ${\phi_A}$ % TFP AR(1)
           phi_R    ${\phi_R}$ % Monetary AR(1)
           phi_Pi ${\phi_{\Pi}}$ ;   % cost push AR(1)


alpha = 0.65;
gamma = 1.5;
beta  = 0.99;
Nbar  = 0.95;    // \overline{N}
kappa = 1.2;
mu    = 6;
rho   = 0.044;

% Additional parameters required for the model but not the steady state

chi   = 0.8;   // Example: for (1 - chi^(-gamma)) in Euler equation
psi   = 60;   // Example: slope parameter in the NK Phillips curve
rho_W = 0.9;   // Example: wage rule inertia
rho_A = 0.5;  // Example: wage rule loading on tfp
delta = 1.5;   // Example: Taylor-type rule coefficient on inflation
phi_A = 0.95;   // Example: AR(1) for TFP
phi_R = 0.95;   // Example: AR(1) for Monetary
phi_Pi = 0.95;   // Example: AR(1) for Cost push



%%%%%%%%%%%% 3. Model Equations (log-linear) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model(linear);


    % 1) M_bar = rho * N_bar
    #Mbar  = rho * Nbar;

    % 2) e_bar = 1 - (1 - rho) * Nbar
    #ebar  = 1 - (1 - rho) * Nbar;

    % 3) v_bar = (Mbar / (ebar^alpha))^(1/(1 - alpha))
    #vbar  = ( Mbar / (ebar^alpha) )^( 1/(1 - alpha) );

    % 4) f_bar = Mbar / vbar
    #fbar  = Mbar / vbar;

    % 5) eta_bar = Mbar / ebar
    #etabar = Mbar / ebar;

    % 6) MC_bar = 1 - 1/mu
    #MCbar = 1 - 1/mu;

    % 7) W_bar = MC_bar - (kappa / fbar) + beta*(1-rho)*(kappa / fbar)
    #Wbar = MCbar - (kappa / fbar) + beta * (1 - rho) * (kappa / fbar);

    % 8) R_bar = (Wbar^(-gamma)) / [ beta * ( (1 - rho*(1 - etabar)) * Wbar^(-gamma)
    %                       + rho*(1 - etabar)*(chi*Wbar)^(-gamma) ) ]

    #Rbar = Wbar^(-gamma) / ( beta * ( (1 - rho*(1 - etabar)) * Wbar^(-gamma)  + rho*(1 - etabar)*(chi*Wbar)^(-gamma) ) );

    % (1) Euler equation 
    -gamma*w = r - pi(+1)  - gamma*beta*Rbar*(1 - rho*(1 - etabar))*w(+1) + beta*Rbar*rho*etabar*(1 - chi^(-gamma))*eta(+1);

    % (2) Production function
    y = a + n;

    % (3) Law of motion for employment
    n = (1 - rho)*n(-1) + rho*m;

    % (4) Job searchers
    e = - ( (1 - rho)*Nbar / ebar ) * n(-1);

    % (5) Matching function
    m = alpha*e + (1 - alpha)*v;

    % (6) Job filling rate
    f = m - v;

    % (7) Job finding rate
    eta = m - e;

    % (8) Marginal cost
    mc = -a + (Wbar/MCbar)*w - (kappa/(fbar*MCbar))*f + beta*(1 - rho)*(kappa/(fbar*MCbar))*f(+1);

    % (9) New Keynesian Phillips curve
    pi = beta*pi(+1) + (mu - 1)/psi * mc + u_Pi;

    % (10) Wage rule
    w = rho_W*w(-1) + (1 - rho_W)*rho_A*a;

    % (11) Monetary policy rule
    r = delta*pi + u_R;

    % (12) TFP process
    a = phi_A*a(-1) + eA - eA4(-4) - eA8(-8);

    % (13) Monetary process
    u_R = phi_R*u_R(-1) + eR ;

    % (14) Cost push process
    u_Pi = phi_Pi*u_Pi(-1) + ePi ;

    % (15) Output growth
    dy = y-y(-1);

    % (16) Unemployment
    u = -n;

end; 

%%%%%%%%%%%% 4. Initial Values for the Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initval;
  w   = 0;  pi  = 0;  r   = 0;  y   = 0;
  n   = 0;  e   = 0;  v   = 0;  m   = 0;
  f   = 0;  eta = 0;  mc  = 0;  a   = 0;
  dy = 0; u = 0;
  eA  = 0;  eA4 = 0; eA8=0;
end;

steady;  % Solve for the steady state
check;   % Check for eigenvalues, etc.

stoch_simul(irf=20);

%%%%%%%%%%%% 5. Declaring Observables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We assume your data are already in the correct form (demeaned or in log-deviations)
% so that dy, r, pi, u match the model variables directly.

varobs dy r pi n;

%%%%%%%%%%%% 6. Estimation: Specifying Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
  % Parameter  , initial guess ,  prior dist  ,  prior mean , prior std
    alpha       ,    0.65       ,  beta_pdf    ,   0.50      , 0.10;
   gamma        ,    2.00       ,  gamma_pdf  ,   2.00      , 0.50;
   psi         ,    60.0       ,  gamma_pdf   ,   75.0      , 25.0;
    kappa       ,    1.2        ,  gamma_pdf   ,   1.0       , 0.3;
   delta       ,    1.5        ,  normal_pdf  ,   1.5       , 0.2;
  % chi         ,    0.8       ,  beta_pdf    ,   0.85      , 0.10;
  % rho         ,    0.04       ,  beta_pdf    ,   0.04      , 0.01;
  phi_A       ,    0.95       ,  beta_pdf    ,   0.50      , 0.20;
  phi_Pi       ,   0.95       ,  beta_pdf    ,   0.50      , 0.20;
  phi_R       ,    0.3       ,  beta_pdf    ,   0.50      , 0.20;
  rho_W       ,    0.3       ,  beta_pdf    ,   0.50      , 0.20;
  rho_A       ,    0.3       ,  normal_pdf    ,   0.50      , 0.20;


  % Estimate the standard errors of the shocks 
  %  (We assume initial guesses ~0.01, and for example an inverse gamma prior.)
  stderr eA   ,  0.01 , inv_gamma_pdf ,  0.001 , 0.1;
  stderr eA4  ,  0.01 , inv_gamma_pdf ,  0.001 , 0.1;
  stderr eA8  ,  0.01 , inv_gamma_pdf ,  0.001 , 0.1;
  stderr ePi  ,  0.01 , inv_gamma_pdf ,  0.001 , 0.1;
  stderr eR   ,  0.01 , inv_gamma_pdf ,  0.001 , 0.1;
end;

identification;

%%%%%%%%%%%% 7. Estimation Command %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimation( datafile  = NSPSLM_estimation_data,  % or .mat, .csv, etc.
            first_obs = 1,
             mode_compute = 6,       % MCMC or optimizer (6=csminwel, or 4=FminSearchSimplex, etc.)
             mh_tune_jscale,
            % mh_replic = 0,          % 0 => do a mode-finder only => ML or Laplace approx
            % Uncomment next line to do full Bayesian with MCMC:
             mh_replic = 100000, 
             mh_nblocks=1, 
             mh_drop=0.2,
             mode_check,  % helps debug the posterior shape
            graph_format = (pdf),
            plot_priors = 1 ,
            bayesian_irf,
            tex,
                 irf=20) a y r pi u v;

write_latex_prior_table;  

shock_decomposition dy n r pi a;

% plot_shock_decomposition y n r pi a;

stoch_simul(irf=20) y n r pi a;

horizon = 1:options_.irf;

% Create a new figure with subplots
figure('Name','IRFs to eA4 Shock','NumberTitle','off');

%--------------------------------------------------------------------------
% 1) TFP (a)
%--------------------------------------------------------------------------
subplot(2,3,1)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.a_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.a_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.a_eA4, 'r--', 'LineWidth', 2); 
title('TFP','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 2) Output (y)
%--------------------------------------------------------------------------
subplot(2,3,2)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.y_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.y_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.y_eA4, 'r--', 'LineWidth', 2); 
title('Output','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 3) Unemployment (u)
%--------------------------------------------------------------------------
subplot(2,3,3)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.u_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.u_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.u_eA4, 'r--', 'LineWidth', 2); 
title('Unemployment','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 4) Vacancies (v)
%--------------------------------------------------------------------------
subplot(2,3,4)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.v_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.v_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.v_eA4, 'r--', 'LineWidth', 2); 
title('Vacancies','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 5) Nominal Interest Rate (r)
%--------------------------------------------------------------------------
subplot(2,3,5)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.r_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.r_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.r_eA4, 'r--', 'LineWidth', 2); 
title('Interest Rate','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 6) Inflation (pi)
%--------------------------------------------------------------------------
subplot(2,3,6)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.pi_eA4, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.pi_eA4, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.pi_eA4, 'r--', 'LineWidth', 2); 
title('Inflation','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

% Adjust layout
sgtitle('Estimated IRFs to the One-Year News Shock','Interpreter','latex','FontSize',14);

print('Figure11A','-dpng');

% Create a new figure with subplots
figure('Name','IRFs to eA4 Shock','NumberTitle','off');

%--------------------------------------------------------------------------
% 1) TFP (a)
%--------------------------------------------------------------------------
subplot(2,3,1)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.a_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.a_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.a_eA8, 'r--', 'LineWidth', 2); 
title('TFP','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 2) Output (y)
%--------------------------------------------------------------------------
subplot(2,3,2)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.y_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.y_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.y_eA8, 'r--', 'LineWidth', 2); 
title('Output','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 3) Unemployment (u)
%--------------------------------------------------------------------------
subplot(2,3,3)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.u_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.u_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.u_eA8, 'r--', 'LineWidth', 2); 
title('Unemployment','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 4) Vacancies (v)
%--------------------------------------------------------------------------
subplot(2,3,4)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.v_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.v_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.v_eA8, 'r--', 'LineWidth', 2); 
title('Vacancies','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
ylim([1.2*min(100.*oo_.PosteriorIRF.dsge.Median.v_eA8) 0]);


%--------------------------------------------------------------------------
% 5) Nominal Interest Rate (r)
%--------------------------------------------------------------------------
subplot(2,3,5)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.r_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.r_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.r_eA8, 'r--', 'LineWidth', 2); 
title('Interest Rate','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

%--------------------------------------------------------------------------
% 6) Inflation (pi)
%--------------------------------------------------------------------------
subplot(2,3,6)
plot(horizon, 100.*oo_.PosteriorIRF.dsge.Median.pi_eA8, 'r-', 'LineWidth', 2); 
hold on 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDinf.pi_eA8, 'r--', 'LineWidth', 2); 
plot(horizon, 100.*oo_.PosteriorIRF.dsge.HPDsup.pi_eA8, 'r--', 'LineWidth', 2); 
title('Inflation','Interpreter','latex','FontSize',12);
xlabel('Time','Interpreter','latex','FontSize',12);
ylabel('Deviation from SS','Interpreter','latex','FontSize',12);
grid on;
yl = ylim;
if yl(1) > 0
    ylim([0 yl(2)]);
end

% Adjust layout
sgtitle('Estimated IRFs to the Two-Year News Shock','Interpreter','latex','FontSize',14);

print('Figure11B','-dpng');


