function [F,jac]=Funz_jacob_plasma(x, x0,v_bc,n_bc,p_bc,dt,r,mun,mup,eps,alpha,S,Vth,gen)
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

Ar_full = plasma_rhs3(r, [v_bc(1); v; v_bc(2)], mun, alpha, Vth, -1);

% R_full = Plasma_rhs(r, [v_bc(1); v; v_bc(2)], [n_bc(1); n; n_bc(2)], mun, alpha, Vth, -1);

% NB Se uso invece Plasma_rhs1 viene restituito direttamente il vettore
% moltiplicato per n, al quale devo togliere il primo e l'ultimo termine,
% ma è già comprensivo della correzione con i termini di bordo, dunque devo
% levarli da sotto (riga 66)


A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);
Ar = Ar_full(2:end-1,2:end-1);
% R = R_full(2:end-1);


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
% int_S=[((r(1:end-1)+r(2:end))/2).^2;0]-[0; ((r(1:end-1)+r(2:end))/2).^2];
% S = S.*int_S(2:end-1);

% Includes the S term (for when using non constant generation term)
% rhs = [zeros(lr-2,1);  dt*M*S + M*n0  ; dt*M*S + M*p0];

% doesn't include S (for when using constant integration term)
rhs = [zeros(lr-2,1);  M*n0  ; M*p0];


bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];
full_rhs = rhs - bounds;


% for plasma_rhs1
% F = NL*x + [zeros(lr-2,1); (-dt*R); (-dt*R)] - full_rhs ;

% for plasma_rhs2 and 3 (best working so far?)
% F = NL*x + [zeros(lr-2,1); -dt*(Ar*n+Ar_bc*n_bc); -dt*(Ar*n+Ar_bc*n_bc)] - full_rhs ;

% constant throughout (useful for testing code only, very unphysical)
% F = NL*x + [zeros(lr-2,1); (-dt*M*ones(lr-2,1)); (-dt*M*ones(lr-2,1))] - full_rhs ;

% with gen term defined in dati_plasma
F = NL*x + [zeros(lr-2,1); (-dt*M*gen); (-dt*M*gen)] - full_rhs ;



if nargout>1
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
