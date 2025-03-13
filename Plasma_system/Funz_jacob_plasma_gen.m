function [F,jac]=Funz_jacob_plasma_gen(x, x0,v_bc,n_bc,p_bc,dt,r,mun,mup,eps,alpha,S,Vth,gen)
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


%% BC CORRECTION ----------------------------------------------------------
A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);


% getting the first and last elements of the matrices 
A_bc = A_full(2:end-1,[1 end]);
An_bc = An_full(2:end-1,[1 end]);
Ap_bc = Ap_full(2:end-1,[1 end]);


%% FULL MATRIX CONSTRUCTION -----------------------------------------------
zeri = zeros(lr-2);

NL = [A M -M; 
    zeri (M+dt*An) zeri; 
    zeri zeri (M+dt*Ap)];

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);


%% RHS AND BOUNDS ---------------------------------------------------------
% doesn't include S (for when using constant integration term)
rhs = [zeros(lr-2,1);  M*n0  ; M*p0];

bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];
full_rhs = rhs - bounds;


%% FULL SYSTEM ------------------------------------------------------------
% with gen term defined in dati_plasma
F = NL*x + [zeros(lr-2,1); (-dt*M*gen); (-dt*M*gen)] - full_rhs ;


%% JACOBIAN ---------------------------------------------------------------
if nargout>1
    vJ=[v_bc(1); v; v_bc(2)];    % Li sovrascrivo perché non servono più
    nJ=[n_bc(1); n; n_bc(2)];
    pJ=[p_bc(1); p; p_bc(2)];

    dAn = zeros(lr, lr);    % NB Si calcola anche la derivata in v0 e v(lr)( condiz al bordo), ma tanto sotto le rimuovo
    dAp = zeros(lr, lr);    

    for i=2:lr-1

        log_pos = 1/log(r(i)/r(i+1));
        log_neg = 1/log(r(i-1)/r(i));
        
        up_neg = (vJ(i+1) - vJ(i))/Vth;
        up_pos = (vJ(i) - vJ(i+1))/Vth;
        dw_neg = (vJ(i) - vJ(i-1))/Vth;
        dw_pos = (vJ(i-1) - vJ(i))/Vth;


        dAn(i,i-1) = log_neg * (nJ(i)*DB(dw_neg) + nJ(i-1)*DB(dw_pos));
        dAn(i,i) = log_pos * (-nJ(i+1)*DB(up_neg) - nJ(i)*DB(up_pos)) +...
                        log_neg * (-nJ(i)*DB(dw_neg) - nJ(i-1)*DB(dw_pos));
        dAn(i,i+1) = log_pos * (nJ(i+1)*DB(up_neg) + nJ(i)*DB(up_pos));


        dAp(i,i-1) = log_neg * (-pJ(i)*DB(-dw_neg) - pJ(i-1)*DB(-dw_pos));
        dAp(i,i) = log_pos * (pJ(i+1)*DB(-up_neg) + pJ(i)*DB(-up_pos)) +...
                        log_neg * (pJ(i)*DB(-dw_neg) + pJ(i-1)*DB(-dw_pos));
        dAp(i,i+1) = log_pos * (-pJ(i+1) * DB(-up_neg) - pJ(i)*DB(-up_pos));

    end
    
    % potrebbero esserci dei termini che dividono dAn e dAp dai vol finiti
    
    J11 = A;
    J12 = M;
    J13 = -M;
    J21 = dt*mun*dAn(2:end-1, 2:end-1);
    J22 = M + dt*An;
    J23 = zeri;
    J31 = dt*mup*dAp(2:end-1, 2:end-1); 
    J32 = zeri;
    J33 = M + dt*Ap;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];

end


end
