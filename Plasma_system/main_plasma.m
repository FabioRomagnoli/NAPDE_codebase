close all;
clear;

addpath('.\Plasma_system\');
load('.\Plasma_system\ExpData.mat')

idx = 14;   % 14
Iz = xxZheng(idx,"I").(1) * 10^(-4);     % microV/cm converted to V/m
Vz = xxZheng(idx,"V").(1) * 10^(3);      % kV converted to V

% Hyper parameters
K = 20;                      % Time grid points
lr = 101;
dt = 1e-4;                    % Time separation  [s]
Vsrt = Vz;                   % Voltage at r=1 and t=1  [V]
Vend = Vz;                   % Ending voltage at r=1 and t=K*dt  [V]
S = 1e9;                      % random constant  [?]
N = 1e7;                      % density constant [m-3]
alphaZ = 1.13968e+04;  % 54;



Dati_plasma;


%% SOLVE ------------------------------------------------------------------
% true for constant generation term for alpha*Jn gen term
solveConstGen = true;

% load steady state solution of const gen term
% if ~solveConstGen
%     load('.\Plasma_system\Xsol.mat');
%     X(:,1) = Xsol;
% end

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
Ar_full = Plasma_rhs2(r, vEnd, mun, 1, Vth, -1);
reactionTerm = Ar_full*nEnd;

% Integration of J with Plasma_rhs2 
intJAll = sum(reactionTerm);
intJ1 = sum(reactionTerm(1:idx));
% intJ2 = sum(reactionTerm(idx+1:end));

% Integral of const generation term with mass matrix 
intGen = sum(M_full*genfull);

% Alpha is ratio of the two integrals over the ionization areas 
alphaAll = intGen/intJAll;
alphaSG = intGen/intJ1; % senza gobba

fprintf('alphaAll = %.3f\n', alphaAll);
fprintf('alphaSG =  %.3f\n', alphaSG);
fprintf('alphaDif = %.3f\n', abs(alphaAll - alphaSG));



%% PLOTTING GEN TERM ------------------------------------------------------
M_full = ax_mass(r, 1);

% Test for integration over i-1/2 to i+1/2 with Mass matrix
x_medi = (r(1:end-1)+r(2:end))/2;
dd = [0; diff(x_medi)/2] + [diff(x_medi)/2; 0];


alphaPlot = alphaZ;    
R_full = Plasma_rhs(r, vEnd, nEnd, mun, alphaPlot, Vth, -1);
Ar_full = Plasma_rhs2(r, vEnd, mun, alphaPlot, Vth, -1);
Jn = dd.*Comp_current(r,mun,q,vEnd,Vth,-1,nEnd);

figure()
title('Comparison of Integral of different generation terms')
hold on; 
plot(x_medi,alphaPlot*Jn/(2*pi*q),"b-o", 'DisplayName', 'Jn');
plot(r, M_full*genfull, "k-s", 'DisplayName', 'Const Gen');
plot(r, Ar_full*nEnd,"-*",'DisplayName', 'Rhs2' );    
plot(r, R_full,"r-x", 'DisplayName', 'Rhs1'); 
set(gca, 'YScale', 'log') % Change y-axis to log scale
set(gca, 'XScale', 'log') % Change y-axis to log scale
legend('Location', 'best'); 
hold off;
grid on;



%% DIMENSIONAL ANALYSIS

% Av_grad = ax_gradient(r); % Computes the gradient with FD of second order
% JnTerm1 = Av_grad*nEnd;   % dn/dx
% JnTerm2 = -(nEnd.*Av_grad*vEnd)/Vth; % -n*dv/dx
% JnTest = q*Vth*mun*abs(JnTerm1 - JnTerm2);
% 
% figure;
% plot(r,r.*JnTest,"k-o");
% set(gca, 'YScale', 'log') % Change y-axis to log scale
% set(gca, 'XScale', 'log') % Change y-axis to log scale
% grid on;

