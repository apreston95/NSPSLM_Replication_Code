// Three-equation New Keynesian model with TFP news shock (Dynare 6.4)

// -----------------------------------------------------------------------------
// Declarations
// -----------------------------------------------------------------------------
var y y_nat y_gap pi i r_n a n w;
varexo eps_news;

parameters beta sigma psi kappa phi_pi phi_y rho_a sigma_a sigma_news;

// -----------------------------------------------------------------------------
// Calibration
// -----------------------------------------------------------------------------
beta       = 0.99;
sigma      = 2;       // Inverse IES
psi        = 1;       // Inverse Frisch elasticity
kappa      = 0.085;
phi_pi     = 1.5;
phi_y      = 0;
rho_a      = 0.99;
sigma_a    = 0.01;    // not used directly (no contemporaneous eps_a here)
sigma_news = 0.01;

// -----------------------------------------------------------------------------
// Model
// -----------------------------------------------------------------------------
model;
// Production (all in logs/deviations)
y = a + n;

// Natural output
y_nat = ((1+psi)/(sigma+psi)) * a;

// Output gap
y_gap = y - y_nat;

// Labour supply (in logs)
w = sigma*y + psi*n;

// IS curve
y = y(+1) - (1/sigma) * ( i - pi(+1) - r_n );

// NK Phillips curve
pi = beta*pi(+1) + kappa*y_gap;

// Taylor rule
i = phi_pi*pi + phi_y*y_gap;

// Natural real rate
r_n = sigma*( y_nat(+1) - y_nat );

// TFP with 4-quarter-ahead news (realisation occurs four periods after signal)
// Using lagged eps_news means a shock at t-4 moves a at t.
a = rho_a*a(-1) - eps_news(-4);
end;

// -----------------------------------------------------------------------------
// Steady state and initial values (zero steady state in log deviations)
// -----------------------------------------------------------------------------
initval;
y = 0;
y_nat = 0;
y_gap = 0;
pi = 0;
i = 0;
r_n = 0;
a = 0;
n = 0;
w = 0;
end;

steady;
check;

// -----------------------------------------------------------------------------
// Shock variances
// -----------------------------------------------------------------------------
shocks;
var eps_news; stderr sigma_news;
end;

// -----------------------------------------------------------------------------
// Simulation
// -----------------------------------------------------------------------------
stoch_simul(order=1, irf=21) y y_nat y_gap pi i a n w;

// -----------------------------------------------------------------------------
// MATLAB plotting
// -----------------------------------------------------------------------------
% IRF plotting (uses oo_.irfs fields like y_eps_news, i_eps_news, n_eps_news)
figure;

darkGreen = [0, 0.5, 0];
black     = [0, 0, 0];

% Output IRF
subplot(1,3,1);
plot(0:20, oo_.irfs.y_eps_news*100, 'LineWidth', 2, 'Color', darkGreen);
title('Output', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\%$ Deviation', 'Interpreter', 'latex');
yline(0, 'LineWidth', 2, 'Color', black);
ylim([-3, 0.5]);

% Nominal interest rate IRF
subplot(1,3,2);
plot(0:20, oo_.irfs.i_eps_news*100, 'LineWidth', 2, 'Color', darkGreen);
title('Nominal Interest Rate', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\%$ Deviation', 'Interpreter', 'latex');
yline(0, 'LineWidth', 2, 'Color', black);
ylim([-2, 0.5]);

% Employment IRF
subplot(1,3,3);
plot(0:20, oo_.irfs.n_eps_news*100, 'LineWidth', 2, 'Color', darkGreen);
title('Employment', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Quarters', 'Interpreter', 'latex');
ylabel('$\%$ Deviation', 'Interpreter', 'latex');
yline(0, 'LineWidth', 2, 'Color', black);
ylim([-1.5, 0.5]);

print('FigureA9','-dpng');
