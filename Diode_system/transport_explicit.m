function [F, jac]=transport_explicit(x, x0,v_bc,n_bc,p_bc,dt,r,N,mu,eps,ni,Vth,tau)
% v_bc,n_bc e l'altro sono vettori colonna di 2 elementi che contengono le
% condizioni al bordo di quelle var
% NB x va passato colonna, è ridotto, dunque ha lunghezza lr-6
%NB x0 è x_prec
% r è vettore colonna

lr=length(r);
v0 = x0(1:lr-2);
n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);
N=N(2:end-1);

A_full = ax_laplacian (r,eps);
M_full = ax_mass(r, 1);
An_full = ax_dd(r, [v_bc(1); v0; v_bc(2)], mu, Vth, -1);
Ap_full = ax_dd(r, [v_bc(1); v0; v_bc(2)], mu, Vth, 1);% Così il numero di valenza è corretto

A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);
An = An_full(2:end-1,2:end-1);
Ap = Ap_full(2:end-1,2:end-1);

zeri = zeros(lr-2);


NL = [A M -M; 
    zeri M zeri; 
    zeri zeri M];


A_bc = A_full(2:end-1,[1 end]);
An_bc = An_full(2:end-1,[1 end]);
Ap_bc = Ap_full(2:end-1,[1 end]);

rhs = [M*N; M*n0; M*p0];
bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];

full_rhs = rhs - bounds;


F = NL*x - full_rhs + [zeros(lr-2,1); dt*An*n0; dt*Ap*p0] ;


if nargout>1    
    J11 = A;
    J12 = M;
    J13 = -M;
    J21 = zeri;
    J22 = M ;
    J23 = zeri;
    J31 = zeri; 
    J32 = zeri;
    J33 = M;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];
end

end
