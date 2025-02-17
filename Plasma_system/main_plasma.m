close all;
clear;

% Hyper parameters
K = 5;                       % Time grid points
lr = 101;
dt = 1e-4;                    % Time separation  [s]
Vsrt = 1.5e4;                 % Voltage at r=1 and t=1  [V]
Vend = 1.5e4;                 % Ending voltage at r=1 and t=K*dt  [V]
S = 1e9;                      % random constant  [no clue]
N = 1e7;                      % density constant [m-3]

% The voltage at r=end is kept at 0

Dati_plasma;
solve_plasma;
plot_plasma;