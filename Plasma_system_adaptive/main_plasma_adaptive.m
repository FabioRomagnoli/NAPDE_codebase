close all;
clear;
clc;

addpath('.\Plasma_system\');
load('.\Plasma_system\ExpData.mat')

idx = 14;   % 14
Iz = xxZheng(idx,"I").(1) * 10^(-4);     % microV/cm converted to V/m
Vz = xxZheng(idx,"V").(1) * 10^(3);      % kV converted to V

% Hyper parameters
K = 1000001;                     % Time save points
T = 100;                        % Final time [s]
lr = 101;
dt0 = 1e-4;                    % Time separation  [s]
Vsrt = 0;                   % Voltage at r=1 and t=1  [V]
Vend = Vz;                   % Ending voltage at r=1 and t=K*dt  [V]
S = 0;                      % random constant  [?]
N = 1e7;                      % density constant [m-3]
alphaZ = 11586.738;  
tol = 1e-3;


Dati_plasma_adaptive;


%% SOLVE ------------------------------------------------------------------

solve_plasma_adaptive;


%% PLOT -------------------------------------------------------------------
% set to 0 for no plot, to 1 for animation, to K for last plot. 
concentration_plot = 0;
potential_plot = 1;
current_plot = 1;

plot_plasma;

%% POST PROCESSING --------------------------------------------------------
clc

Jn = Comp_current(r,mun,q,vEnd,Vth,-1,nEnd);
Jp = Comp_current(r,mup,q,vEnd,Vth, 1,pEnd);
JJ = Jn + Jp;

% Check if JJ is constant
if std(JJ) / mean(JJ) < 1e-2 
    Ic = mean(JJ);
    Id = abs(Iz - Ic);
    % Display results
    fprintf('Vz = %.4s\n', Vz)
    fprintf('Iz = %.5s\n', Iz);
    fprintf('Ic = %.5s\n', Ic);
    fprintf('Id = %.5s\n', Id);
else
    fprintf('JJ is not constant. \n');
end