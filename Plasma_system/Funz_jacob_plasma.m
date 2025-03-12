function [F,jac]=Funz_jacob_plasma(x, x0,v_bc,n_bc,p_bc,dt,r,mun,mup,eps,alpha,S,Vth, beta, Ei)
% v_bc,n_bc e l'altro sono vettori colonna di 2 elementi che contengono le
% condizioni al bordo di quelle var
% NB x va passato colonna, è ridotto, dunque ha lunghezza lr-6
% NB x0 è x_prec
% r è vettore colonna

%% splitting vector in its components
lr = length(r);
v = x(1:lr-2);
n = x(lr-1:2*lr-4);
p = x(2*lr-3:end);


%% MATRIX DEFINITIONS -----------------------------------------------------
A_full = ax_laplacian (r,eps);
M_full = ax_mass(r, 1);
An_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mun, Vth, -1);
Ap_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mup, Vth, 1); % Così il numero di valenza è corretto

% Av_grad = ax_gradient(r);
% E = -Av_grad*[v_bc(1); v; v_bc(2)];
% alpha = beta .* exp(-Ei ./ E); % Townsend formula

% % Plot results
% figure;
% semilogy(r, alpha, 'b', 'LineWidth', 2); % Log-scale for better visualization
% xlabel('Radius r (m)');
% ylabel('Ionization Coefficient \alpha (1/m)');
% title('Townsend Ionization Coefficient \alpha(r)');
% grid on;

%% GENERATION TERM --------------------------------------------------------
Ar_full = Plasma_rhs2(r, [v_bc(1); v; v_bc(2)], mun, alpha, Vth, -1);


%% BC CORRECTION ----------------------------------------------------------
A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);
Ar = Ar_full(2:end-1,2:end-1);


% getting the first and last elements of the matrices 
A_bc = A_full(2:end-1,[1 end]);
An_bc = An_full(2:end-1,[1 end]);
Ap_bc = Ap_full(2:end-1,[1 end]);
Ar_bc = Ar_full(2:end-1,[1 end]);


%% FULL MATRIX CONSTRUCTION -----------------------------------------------
zeri = zeros(lr-2);

NL = [A M -M; 
    zeri (M+dt*An) zeri; 
    zeri zeri (M+dt*Ap)];

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);


%% RHS AND BOUNDS ---------------------------------------------------------
% Includes the S term (for when using non constant generation term)
rhs = [zeros(lr-2,1);  dt*M*S + M*n0  ; dt*M*S + M*p0];

bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];
full_rhs = rhs - bounds;


%% FULL SYSTEM ------------------------------------------------------------

Rhs_red = zeros(lr-2,1);
Rhs_red(1:36) = Ar(1:36,1:36)*n(1:36);

% for plasma_rhs2 and 3 (best working so far?)
F = NL*x + [zeros(lr-2,1); -dt*(Rhs_red+Ar_bc*n_bc); -dt*(Ar*n+Ar_bc*n_bc)] - full_rhs ;


%% JACOBIAN ---------------------------------------------------------------
if nargout>1
    % Non-constant reaction term
    dRn = dt*Ar;  %derivata di R rispetto a n

    % Constant reaction term
    % dRn = zeri;  %derivata di R rispetto a n

    % dRp is zero regardles since the generation term is based only on n

    v=[v_bc(1); v; v_bc(2)];    % Li sovrascrivo perché non servono più
    n=[n_bc(1); n; n_bc(2)];
    p=[p_bc(1); p; p_bc(2)];

    dAn = zeros(lr, lr);    % NB Si calcola anche la derivata in v0 e v(lr)( condiz al bordo), ma tanto sotto le rimuovo
    dAp = zeros(lr, lr);    

    for i=2:lr-1

        log_pos = 1/log(r(i)/r(i+1));
        log_neg = 1/log(r(i-1)/r(i));
        
        up_neg = (v(i+1) - v(i))/Vth;
        up_pos = (v(i) - v(i+1))/Vth;
        dw_neg = (v(i) - v(i-1))/Vth;
        dw_pos = (v(i-1) - v(i))/Vth;


        dAn(i,i-1) = log_neg * (n(i)*DB(dw_neg) + n(i-1)*DB(dw_pos));
        dAn(i,i) = log_pos * (-n(i+1)*DB(up_neg) - n(i)*DB(up_pos)) +...
                        log_neg * (-n(i)*DB(dw_neg) - n(i-1)*DB(dw_pos));
        dAn(i,i+1) = log_pos * (n(i+1)*DB(up_neg) + n(i)*DB(up_pos));


        dAp(i,i-1) = log_neg * (-p(i)*DB(-dw_neg) - p(i-1)*DB(-dw_pos));
        dAp(i,i) = log_pos * (p(i+1)*DB(-up_neg) + p(i)*DB(-up_pos)) +...
                        log_neg * (p(i)*DB(-dw_neg) + p(i-1)*DB(-dw_pos));
        dAp(i,i+1) = log_pos * (-p(i+1) * DB(-up_neg) - p(i)*DB(-up_pos));

    end
    
    % potrebbero esserci dei termini che dividono dAn e dAp dai vol finiti
    
    J11 = A;
    J12 = M;
    J13 = -M;
    J21 = dt*mun*dAn(2:end-1, 2:end-1);
    J22 = M + dt*An - dRn;
    J23 = zeri;
    J31 = dt*mup*dAp(2:end-1, 2:end-1); 
    J32 = -dRn;
    J33 = M + dt*Ap;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];

end


end
