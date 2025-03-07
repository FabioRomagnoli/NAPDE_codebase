close all;
clear;

% Hyper parameters
K = 50;                      % Time grid points
lr = 201;
dt = 1e-4;                    % Time separation  [s]
Vsrt = 3.017e4;               % Voltage at r=1 and t=1  [V]
Vend = 3.017e4;               % Ending voltage at r=1 and t=K*dt  [V]
S = 1e9;                      % random constant  [?]
N = 1e7;                      % density constant [m-3]

% Solver options:
% Nothing provided
% options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');

% Full Jacobian
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter'); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);

% The voltage at r=end is kept at 0
Dati_plasma;
solve_plasma;
plot_plasma;