%% MESH -------------------------------------------------------------------
% Non Uniform mesh
r0 = 700e-6;
r1 = r0 + 10.35e-2;

% Mesh 
% delta   = 2.4577e-05; % original value  
delta   = 2.4577e-07;

alpha_msh  = -0.1;
msh = CreateTanhMsh(lr, r0, r1-r0, delta, alpha_msh);
r = msh.x;

%% Adaptive
tsave = linspace (0, T, K);

%% PARAMETERS -------------------------------------------------------------
epsilon = 8.8e-12;             % Permittivity              [C]/([V][m])
mup = 2e-4;                    % Diffusivity coefficients  [m2]/([s][V])
mun = 1e2 * mup; 

q = 1.6e-19;                   % Charge                    [C] 
Vth = 26e-3;                   % Thermal voltage           [V]

S = S*ones(lr-2,1);   

%% INITIAL CONDITIONS -----------------------------------------------------
n0 = ones(lr,1) * N;
p0 = ones(lr,1) * N;

% solving for v0 inside domain 
F = ax_laplacian (r, epsilon);    % Perché cambia il valore del residuo (quello restituito da fsolve) in base a se metto 1 o epsilon? l'rhs è 0!
v0 = zeros(size(r));
v0(1) = Vsrt;               % can be inbetween 0 < v(0) < 5e4   [V]
v0(end) = 0;
% Siccome p-n+N fa 0 allora al rhs c'è solo la correzione con le condizioni al bordo
v0(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * v0([1 end]) + q * (p0(2:lr-1) - n0(2:lr-1)));

% Border conditions
v_bc = [v0(1); v0(end)]; % Nb v_bc(:,1) sono le cond al bordo al passo temp 2
n_bc = [n0(1); n0(end)];
p_bc = [p0(1); p0(end)];
      

%% ADIMENSIONALIZATION ----------------------------------------------------
% Scaling factors
xbar = r1-r0;                    % Lenght of device
nbar = N;                        % This N doesn't have the same meaning as the diode problem
Vbar = Vend;
mubar= max(mun,mup);
tbar = 1/(mubar*Vbar/xbar^2);  
Jbar = mubar*Vbar*nbar/xbar;

% Scaling procedure
epsin = epsilon*Vbar/(q*nbar*xbar^2);     % adimensionale (squared normalized D. L.)
Vthin = Vth / Vbar;
Tin = T/tbar;
dt0in = dt0/tbar; 
rin = r/xbar;
p0in = p0/nbar;
n0in = n0/nbar;
v0in = v0/Vbar;

Sin = S * xbar^2/(mubar*Vbar*nbar);
alphain = alphaZ * xbar;

v_bcin = v_bc/Vbar;
n_bcin = n_bc/nbar;
p_bcin = p_bc/nbar;

munin=mun/mubar; 
mupin=mup/mubar;

Vsrtin = Vsrt/Vbar;
Vendin = Vend/Vbar;

tsavein = tsave/tbar;

%% VECTOR BUILDING --------------------------------------------------------
x0in = [v0in; n0in; p0in];

X = zeros(3*lr,K);
X(:,1) = x0in;
 
X([1 lr]) = v_bcin;
X([lr+1 2*lr]) = n_bcin;
X([2*lr+1 3*lr]) = p_bcin;

