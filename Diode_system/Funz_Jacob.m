function [F, jac]=Funz_Jacob(x, x0,v_bc,n_bc,p_bc,dt,r,N,mu,eps,ni,Vth,tau)
% v_bc,n_bc e l'altro sono vettori colonna di 2 elementi che contengono le
% condizioni al bordo di quelle var
% NB x va passato colonna, è ridotto, dunque ha lunghezza lr-6
%NB x0 è x_prec
% r è vettore colonna

lr=length(r);
v = x(1:lr-2);
n = x(lr-1:2*lr-4);
p = x(2*lr-3:end);
N=N(2:end-1);

A_full = ax_laplacian (r,eps);
M_full = ax_mass(r, 1);
An_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mu, Vth, -1);
Ap_full = ax_dd(r, [v_bc(1); v; v_bc(2)], mu, Vth, 1);% Così il numero di valenza è corretto

A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);

zeri = zeros(lr-2);


NL = [A M -M; 
    zeri (M+dt*An) zeri; 
    zeri zeri (M+dt*Ap)];

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);

R = @(n,p) (ni^2-n.*p)./(tau*(n+p));


A_bc = A_full(2:end-1,[1 end]);
An_bc = An_full(2:end-1,[1 end]);
Ap_bc = Ap_full(2:end-1,[1 end]);

rhs = [M*N; M*n0; M*p0];
bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];

full_rhs = rhs - bounds;


F = NL*x + [zeros(lr-2,1); (-dt*M*R(n,p)); (-dt*M*R(n,p))] - full_rhs ;


if nargout>1
    dR_dn=@(n,p) -(p.^2 + ni^2)./(tau*((n+p).^2));
    dR_dp=@(n,p) -(n.^2 + ni^2)./(tau*((n+p).^2));

    dRn = diag(dt*M*dR_dn(n,p));  %derivata di R rispetto a n, 
    % NB Va sommata perché ci sarebbe un meno moltiplicato a tutto che non ho messo
    dRp = diag(dt*M*dR_dp(n,p)); 

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
    J21 = dt*mu*dAn(2:end-1, 2:end-1);
    J22 = M + dt*An - dRn;
    J23 = -dRp;
    J31 = dt*mu*dAp(2:end-1, 2:end-1); 
    J32 = -dRn;
    J33 = M + dt*Ap - dRp;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];

end

end
