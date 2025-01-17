r = linspace(1,1.000001,lr)';   % grid points       [m]

tau = 1e-11;                     % Time separation   [s]
T = 10*tau;                     % Time span covered [s]

% tempo di transito per diffusion

ni = 6.14e+15;                    % Intrinsic concentration
% epsilon=8.8e-12;                % Permittivity              [C]/[N][m2]
epsilon=1.035e-10;                % Permittivity              [C]/[N][m2]

mu = 0.1;                       % Diffusivity coefficients  [m2]/[s]
q = 1.6e-19;                    % Charge                    [C] 
Vth = 26e-3;                    % Thermal voltage           [V]
N = N*ones(size(r));        % N density of the gas      [m-3]

% Scaling factors
xbar = 1e-6;                    % Lenght of device
nbar = norm(N,'inf');
Vbar = Vth;
tbar = tau;


%% Initial Conditions
%% Caso 1
% La cond iniziale comprende anche i termini di bordo
% Valore esatto dell'n (si trova dal fatto che p=ni^2/n
% e sostituendo nell'eq di neutralità p-n+N)

% n0 = (N+sqrt(N.^2+ni^2*4))/2;   
% p0 = ni^2./n0; 

%% Caso 2
% Nd  = 1e24;
% Na = 1e22;
% 
% N = Nd * (r <= 1 + 500e-9) -...
%     Na * (r > 1 + 500e-9);
% 
% n0 = (N+sqrt(N.^2+ni^2*4))/2;% Valore esatto dell'n (si trova dal fatto che p=ni^2/n
% p0 = (-N+sqrt(N.^2+4*ni^2))/2;
% p0(N>0) = ni^2./n0(N>0); % Capire perché è necessario, forse è solo un problema di approx numerica, infatti N^2 è molto maggiore di ni^2
% n0(N<0) = ni^2./p0(N<0);


%% Caso 3
Nd  = 1e24;
Na1 = 1e22;
Na2 = 1e24;

N = Nd * (r <= 1+50e-9) ...
    - Na1 * (r > 1+50e-9 & r <=1+ 500e-9) ...
    - Na2 * (r > 1+500e-9);

n0 = (N+sqrt(N.^2+ni^2*4))/2;% Valore esatto dell'n (si trova dal fatto che p=ni^2/n
p0 = (-N+sqrt(N.^2+4*ni^2))/2;
p0(N>0) = ni^2./n0(N>0); % Capire perché è necessario, forse è solo un problema di approx numerica, infatti N^2 è molto maggiore di ni^2
n0(N<0) = ni^2./p0(N<0);


F = ax_laplacian (r, 1);% Perché cambia il valore del residuo (quello restituito da fsolve) in base a se metto 1 o epsilon? l'rhs è 0!
v0 = zeros (size (r));
v0(1) = Vsrt;
v0(end) = Vend;
% Siccome p-n+N fa 0 allora al rhs c'è solo la correzione con le condizioni al bordo
v0(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * v0([1 end]));

x0 = [v0; n0; p0];            


% Border conditions
v_bc = [v0(1)*ones(1,K); v0(end)*ones(1, K)];% Nb v_bc(:,1) sono le cond al bordo al passo temp 2
n_bc = [n0(1)*ones(1,K); n0(end)*ones(1, K)];
p_bc = [p0(1)*ones(1,K); p0(end)*ones(1, K)];


% Scaling procedure
epsin = epsilon*Vbar/(q*nbar*xbar^2);     % adimensionale (squared normalized D. L.)
tauin = tau/tbar;
niin = ni/nbar;
Vthin = Vth / Vbar;
Tin = T/tbar;
dtin = dt/tbar;
xin = r/xbar;
Nin  = N/nbar;
p0in = p0/nbar;
n0in = n0/nbar;
v0in = v0/Vbar;

v_bcin = v_bc/Vbar;
n_bcin = n_bc/nbar;
p_bcin = p_bc/nbar;


x0in = [v0in; n0in; p0in];


muin=mu*tbar*Vbar/xbar^2;

X = zeros(3*lr,K+1);
X(:,1) = x0in;

X([1 lr],2:K+1) = v_bcin;
X([lr+1 2*lr],2:K+1) = n_bcin;
X([2*lr+1 3*lr],2:K+1) = p_bcin;

Jp=Jacob_patt(lr-2);




