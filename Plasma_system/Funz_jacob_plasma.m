function [F]=Funz_jacob_plasma(x, x0,v_bc,n_bc,p_bc,dt,r,mun,mup,eps,alpha,S,Vth)
% v_bc,n_bc e l'altro sono vettori colonna di 2 elementi che contengono le
% condizioni al bordo di quelle var
% NB x va passato colonna, è ridotto, dunque ha lunghezza lr-6
% NB x0 è x_prec
% r è vettore colonna

lr = length(r);
v = x(1:lr-2);
n = x(lr-1:2*lr-4);
p = x(2*lr-3:end);

A_full = ax_laplacian (r,eps);
M_full = ax_mass(r, 1);
An_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mun, Vth, -1);
Ap_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mup, Vth, 1); % Così il numero di valenza è corretto
Ar_full=plasma_rhs3(r, [v_bc(1); v; v_bc(2)], mun, alpha, Vth, -1);
%R_full=Plasma_rhs(r, [v_bc(1); v; v_bc(2)], [n_bc(1); n; n_bc(2)], mun, alpha, Vth, -1);

% NB Se uso invece Plasma_rhs1 viene restituito direttamente il vettore
% moltiplicato per n, al quale devo togliere il primo e l'ultimo termine,
% ma è già comprensivo della correzione con i termini di bordo, dunque devo
% levarli da sotto (riga 66)


A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);
Ar=Ar_full(2:end-1,2:end-1);
%R=R_full(2:end-1);


% 
% R_full = ax_r(r, [v_bc(1); v; v_bc(2)], mun, Vth, -1);
% R=R_full(2:end-1,2:end-1);
% R=alpha*R*n;

zeri = zeros(lr-2);

NL = [A M -M; 
    zeri (M+dt*An) zeri; 
    zeri zeri (M+dt*Ap)];

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);


%% Matrix construction

A_bc = A_full(2:end-1,[1 end]);
An_bc = An_full(2:end-1,[1 end]);
Ap_bc = Ap_full(2:end-1,[1 end]);
Ar_bc = Ar_full(2:end-1,[1 end]);

% Forse S dovrebbe essere integrato così anziché fare dt*M*S
%int_S=[((r(1:end-1)+r(2:end))/2).^2;0]-[0; ((r(1:end-1)+r(2:end))/2).^2];

%S=S.*int_S(2:end-1);


rhs = [zeros(lr-2,1);  dt*M*S + M*n0  ; dt*M*S + M*p0];
% rhs = [zeros(lr-2,1);  M*n0  ; M*p0];

bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];

full_rhs = rhs - bounds;

% Jbar = mubar*Vbar*nbar/xbar;


gen = zeros(lr-2,1);
gen(r*0.1035<=800e-6) = 2.652582384864923e+21/5.601064202198418e+11;

%F = NL*x + [zeros(lr-2,1); (-dt*R); (-dt*R)] - full_rhs ;
% F = NL*x + [zeros(lr-2,1); -dt*M*(Ar*n+Ar_bc*n_bc); -dt*M*(Ar*n+Ar_bc*n_bc)] - full_rhs ;
% F = NL*x + [zeros(lr-2,1); (-dt*M*ones(lr-2,1)); (-dt*M*ones(lr-2,1))] - full_rhs ;

F = NL*x + [zeros(lr-2,1); (-dt*M*gen); (-dt*M*gen)] - full_rhs ;


end
