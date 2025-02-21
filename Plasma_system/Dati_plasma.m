% Non Uniform mesh
r0 = 700e-6;
r1 = r0 + 10.35e-2;

% Mesh 
delta   = 2.4577e-05;
% delta   = 2.4577e-03;

alpha_msh  = -0.1;
msh = CreateTanhMsh(lr, r0, r1-r0, delta, alpha_msh);
r = msh.x;
% 
T = K*dt;                      % seconds [s]

epsilon = 8.8e-12;             % Permittivity              [C]/([V][m])
mup = 2e-4;                    % Diffusivity coefficients  [m2]/([s][V])
mun = 1e2 * mup; 

q = 1.6e-19;                   % Charge                    [C] 
Vth = 26e-3;                   % Thermal voltage           [V]

% Plasma related constants
alpha = 7.2e5;         % impact ionization coefficient                    [m]/([s][V])
% alpha = 1.107778321468997e+25;
% alpha = 1e22;

S = S*ones(lr-2,1);   
%E0 = 1e5;                 %  [V]/[m]
%Ei=2.09e7;


% Initial conditions
n0 = ones(lr,1) * N;
p0 = ones(lr,1) * N;


% Initial Condition
F = ax_laplacian (r, epsilon);    % Perché cambia il valore del residuo (quello restituito da fsolve) in base a se metto 1 o epsilon? l'rhs è 0!
v0 = zeros(size(r));
v0(1) = Vsrt;               % can be inbetween 0 < v(0) < 5e4   [V]
v0(end) = 0;
% Siccome p-n+N fa 0 allora al rhs c'è solo la correzione con le condizioni al bordo
v0(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * v0([1 end]) + q * (p0(2:lr-1) - n0(2:lr-1)));


% Scaling factors
xbar = r1-r0;                    % Lenght of device
nbar = N;                        % This N doesn't have the same meaning as the diode problem
Vbar = Vend;
mubar= max(mun,mup);
tbar = 1/(mubar*Vbar/xbar^2) ;  


% Border conditions
v_bc = [linspace(v0(1),Vend,K); v0(end)*ones(1, K)];% Nb v_bc(:,1) sono le cond al bordo al passo temp 2
n_bc = [n0(1)*ones(1,K); n0(end)*ones(1, K)];
p_bc = [p0(1)*ones(1,K); p0(end)*ones(1, K)];


x0 = [v0; n0; p0];            

% Scaling procedure
epsin = epsilon*Vbar/(q*nbar*xbar^2);     % adimensionale (squared normalized D. L.)
Vthin = Vth / Vbar;
Tin = T/tbar;
dtin = dt/tbar; 
rin = r/xbar;
p0in = p0/nbar;
n0in = n0/nbar;
v0in = v0/Vbar;

Sin = S * xbar^2/(mubar*Vbar*nbar);
% alphain = alpha * xbar;
%E0in = E0 * xbar / Vbar; 

v_bcin = v_bc/Vbar;
n_bcin = n_bc/nbar;
p_bcin = p_bc/nbar;

munin=mun/mubar; 
mupin=mup/mubar;

x0in = [v0in; n0in; p0in];

X = zeros(3*lr,K+1);
X(:,1) = x0in;

X([1 lr],2:K+1) = v_bcin;
X([lr+1 2*lr],2:K+1) = n_bcin;
X([2*lr+1 3*lr],2:K+1) = p_bcin;

%% ----------------------------------------

ionization_length = 1.4371e-04;
% ionization_length = 1e-04;

A = ((r0+ionization_length)^2-r0^2)*pi;
% A = ((r0+ionization_length)^2-r0^2)*pi;

I = 0.218e-3;
% I = 0.2e-3;


% gen = zeros(lr-2,K);
% Glin = linspace(0,G,5);
% for k=1:5
%     gen(r <= r0 + ionization_length,k) = Glin(k)/(Jbar/xbar);
% end

G = I/(A*q);
Jbar = mubar*Vbar*nbar/xbar;
gen = zeros(lr-2,1);
gen(r <= r0 + ionization_length) = G/(Jbar/xbar);

% gen(r*0.1035<=800e-6) = 2.652582384864923e+21/5.601064202198418e+11;

% Ic = 2.5e-4;
% j = Ic/(2*pi*r0);
% alphaC = G/j;


% rate di generazione medio : q*G*A = I 
% densità di corrente all'emettitore : j = I/(2*pi*re)
% coefficiente di generazione per impatto : alpha = G/j

alphaAnalytic = (2*pi*r0)/(A*q);
alphain = alphaAnalytic * xbar;
