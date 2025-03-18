close all;
clear;
clc;

addpath('.\Plasma_system\');
load('.\Plasma_system\ExpData.mat')

idx = 14;   % 14
Iz = xxZheng(idx,"I").(1) * 10^(-4);     % microV/cm converted to V/m
Vz = xxZheng(idx,"V").(1) * 10^(3);      % kV converted to V

% Hyper parameters
K = 100;                      % Time grid points
T = 1000;
lr = 101;
dt = 1e-4;                    % Time separation  [s]
Vsrt = Vz;                   % Voltage at r=1 and t=1  [V]
Vend = Vz;                   % Ending voltage at r=1 and t=K*dt  [V]
S = 0;                      % random constant  [?]
N = 1e7;                      % density constant [m-3]
% alphaZ = 11586.738;  % 54
% alphaZ = 11397.060;   
alphaZ = 10;  % 54

Dati_plasma;


%% SOLVE ------------------------------------------------------------------
% true for constant generation term for alpha*Jn gen term
solveConstGen = false;
genplots = false;

% load steady state solution of const gen term
if solveConstGen
    load('.\Plasma_system\Xsol.mat');
    X(:,1) = Xsol;
end

solve_plasma;

% save steady state solution if solved with const gen term
if solveConstGen
    Xsol = X(:,end);
    save(fullfile(".\Plasma_system\", "Xsol.mat"), 'Xsol');
end

%% PLOT -------------------------------------------------------------------
% set to 0 for no plot, to 1 for animation, to K for last plot. 
concentration_plot = 1;
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


%% INTEGRATION ------------------------------------------------------------
genIntFull = M_full*genfull;

R_full = Plasma_rhs(r, vEnd, nEnd, mun, 1, Vth, -1);

% Integration of J with Plasma_rhs2 
intJAll = sum(R_full);
intJ1 = sum(R_full(2:idx));
% intJ2 = sum(reactionTerm(idx+1:end));

% Integral of const generation term with mass matrix 
intGen = sum(genIntFull);

% Alpha is ratio of the two integrals over the ionization areas 
alphaAll = intGen/intJAll;
alphaSG = intGen/intJ1; % senza gobba

fprintf('alphaAll = %.3f\n', alphaAll);
fprintf('alphaSG =  %.3f\n', alphaSG);
fprintf('alphaDif = %.3f\n', abs(alphaAll - alphaSG));


R_full = Plasma_rhs(r, vEnd, nEnd, mun, alphaSG, Vth, -1);
summDiff = sum(genIntFull(2:idx) - R_full(2:idx));
fprintf('summDiff = %.3f\n', summDiff);


figure;
plot(r(2:idx),genIntFull(2:idx) - R_full(2:idx),"k-o");
% set(gca, 'YScale', 'log') % Change y-axis to log scale
set(gca, 'XScale', 'log') % Change y-axis to log scale
grid on;



%% PLOTTING GEN TERM ------------------------------------------------------
% Test for integration over i-1/2 to i+1/2 with Mass matrix
x_medi = (r(1:end-1)+r(2:end))/2;
dd = [0; diff(x_medi)/2] + [diff(x_medi)/2; 0];


alphaPlot = alphaSG;    
R_full = Plasma_rhs(r, vEnd, nEnd, mun, alphaPlot, Vth, -1);
Ar_full = Plasma_rhs2(r, vEnd, mun, alphaPlot, Vth, -1);
Jn = Comp_current(r,mun,q,vEnd,Vth,-1,nEnd);

figure()
title('Comparison of Integral of different generation terms')
hold on; 
% plot(x_medi,alphaPlot*Jn/(2*pi*q),"b-o", 'DisplayName', 'Jn [1/m^2*s]');
% plot(r,genfull,"b-o", 'DisplayName');
plot(r(1:idx), genIntFull(1:idx), "k-s", 'DisplayName', 'Const Gen');
% plot(r, Ar_full*nEnd,"-*",'DisplayName', 'Rhs2' );    
plot(r(1:idx), R_full(1:idx),"r-x", 'DisplayName', 'Rhs1'); 
set(gca, 'YScale', 'log') % Change y-axis to log scale
set(gca, 'XScale', 'log') % Change y-axis to log scale
legend('Location', 'best'); 
ylabel("RHS [1/(ms)]")
hold off;
grid on;



