addpath('.\Plasma_system\');

close all;
clear;

load('.\Plasma_system\ExpData.mat')
load('.\Plasma_system\XSol.mat');

idx = 14;   % 14
Iz = xxZheng(idx,"I").(1) * 10^(-4);     % microV/cm converted to V/m
Vz = xxZheng(idx,"V").(1) * 10^(3);      % kV converted to V

% Hyper parameters
K = 100;                      % Time grid points
lr = 101;
dt = 1e-4;                    % Time separation  [s]
Vsrt = Vz;                   % Voltage at r=1 and t=1  [V]
Vend = Vz;                   % Ending voltage at r=1 and t=K*dt  [V]
S = 0;                      % random constant  [?]
N = 1e7;                      % density constant [m-3]

Dati_plasma;

%% SOLVE ------------------------------------------------------------------
solve_plasma;

%% PLOT -------------------------------------------------------------------
concentration_plot = 1;
potential_plot = 1;
current_plot = 1;

plot_plasma;

%% POST PROCESSING --------------------------------------------------------
clc

% Check if JJ is constant
if std(JJ) / mean(JJ) < 1e-2 
    Ic = mean(JJ);
    Id = abs(Iz - Ic);
    % Display results
    fprintf('Vz = %.2f\n', Vz)
    fprintf('Iz = %.7f\n', Iz);
    fprintf('Ic = %.7f\n', Ic);
    fprintf('Id = %.7f\n', Id);
else
    fprintf('JJ is not constant. \n');
end

%% CALCULATING ALPHA ------------------------------------------------------
M_full = ax_mass(r, 1);
x_medi = (r(1:end-1)+r(2:end))/2;
M_full_centers = ax_mass(x_medi,1);

alphaPlot = 1.1397e+04;

R_full = Plasma_rhs(r, vEnd, nEnd, mun, alphaPlot, Vth, -1);
Ar_full = Plasma_rhs2(r, vEnd, mun, alphaPlot, Vth, -1);
Jn = Comp_current(r,mun,q,vEnd,Vth,-1,nEnd);

figure()
title('Comparison of Integral of different generation terms')
hold on; 
% plot(x_medi,M_full_centers*Jn/(2*pi*q),"b-o", 'DisplayName', 'Jn');
plot(r, M_full*genfull, "k-s", 'DisplayName', 'Const Gen');
plot(r,Ar_full*nEnd,"-*",'DisplayName', 'Rhs2' );    
% plot(r,-R_full,"r-x", 'DisplayName', 'Rhs1'); 
set(gca, 'YScale', 'log') % Change y-axis to log scale
set(gca, 'XScale', 'log') % Change y-axis to log scale
legend('Location', 'best'); 
hold off;
grid on;



%% INTEGRATION ------------------------------------------------------------
reactionTerm = Ar_full*nEnd;

intJAll = sum(reactionTerm);
intJ1 = sum(reactionTerm(1:idx));
% intJ2 = sum(reactionTerm(idx+1:end));

intGen = sum(M_full*genfull);


alphaAll = intGen/intJAll;
alphaSG = intGen/intJ1;



%% DIMENSIOANL ANALYSIS

Av_grad = ax_gradient(r);
JnTerm1 = Av_grad*nEnd;
JnTerm2 = -(nEnd.*Av_grad*vEnd)/Vth;
JnTest = q*Vth*mun*abs(JnTerm1 - JnTerm2);
% JnTest = (nEnd.*Av_grad*vEnd);

% JnTestAbs = (JnTest + abs(JnTest))/2;
figure;
plot(r,r.*JnTest,"k-o");
set(gca, 'YScale', 'log') % Change y-axis to log scale
set(gca, 'XScale', 'log') % Change y-axis to log scale
grid on;

