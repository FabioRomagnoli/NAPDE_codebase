function [F, jac]= reaction_fully_explicit(x, x0,v_bc,n_bc,p_bc,dt,r,N,mu,eps,ni,Vth,tau)
% v_bc,n_bc e l'altro sono vettori colonna di 2 elementi che contengono le
% condizioni al bordo di quelle var
% NB x va passato colonna, è ridotto, dunque ha lunghezza lr-6
%NB x0 è x_prec
% r è vettore colonna

lr=length(r);
v = x(1:lr-2);
n = x(lr-1:2*lr-4);
p = x(2*lr-3:end);
N = N(2:end-1);

A_full = ax_laplacian (r,eps);
M_full = ax_mass(r, 1);

A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);

zeri = zeros(lr-2);

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);

R = (ni^2-n0.*p0)./(tau*(n0+p0));


% NL = [A   zeri  zeri; 
%     zeri  M     zeri; 
%     zeri  zeri  M];

NL = [zeri  zeri  zeri; 
    zeri  M     zeri; 
    zeri  zeri  M];

A_bc = A_full(2:end-1,[1 end]);

% rhs = [M*N-A_bc*v_bc; M*n0+dt*M*R; M*p0+dt*M*R];

rhs = [zeros(lr-2,1); M*n0+dt*M*R; M*p0+dt*M*R];

F = NL*x - rhs + [zeros(lr-2,1); M*n0; -M*p0];

if nargout>1
    J11 = zeri;
    J12 = zeri;
    J13 = zeri;
    J21 = zeri;
    J22 = M;
    J23 = zeri;
    J31 = zeri; 
    J32 = zeri;
    J33 = M;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];

end

end
