clear;

% Hyper parameters
K = 50;                       %Time grid points
lr = 501;
dt = 1e-10;                   % Time separation  [s]
Vsrt = 1.4;                   % Starting voltage (r=1) [V]
Vend = 0;                     % Ending voltage (r=end) [V]
N = 10^22;                    % density constant [m-3]


% Solved with fsolve
Dati_diode;
solve_diode;
plot_diode;

% Solved with custom implementation of Newton
Dati_diode;
solve_diode_newton;
plot_diode;

% Solved with operator splitting
Dati_diode;
solve_diode_split
plot_diode;


