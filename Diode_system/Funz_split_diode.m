function [F, jac]=Funz_split_diode(x, x0,v_bc,n_bc,p_bc,dt,r,N,mu,eps,ni,Vth,tau)
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

A = A_full(2:end-1,2:end-1);
M = M_full(2:end-1,2:end-1);

zeri = zeros(lr-2);

NL = [A M -M; 
    zeri M zeri; 
    zeri zeri M];

n0 = x0(lr-1:2*lr-4);
p0 = x0(2*lr-3:end);

R = @(n,p) (ni^2-n.*p)./(tau*(n+p));

MR = [zeros(lr-2,1); (-dt*M*R(n,p)); (-dt*M*R(n,p))];


A_bc = A_full(2:end-1,[1 end]);


rhs = [M*N; M*n0; M*p0];
bounds = [A_bc*v_bc; zeros(lr-2,1) ; zeros(lr-2,1)];

full_rhs = rhs - bounds;


F = NL*x +  MR - full_rhs ;


if nargout>1
    dR_dn=@(n,p) -(p.^2 + ni^2)./(tau*((n+p).^2));
    dR_dp=@(n,p) -(n.^2 + ni^2)./(tau*((n+p).^2));

    dRn = diag(dt*M*dR_dn(n,p));  %derivata di R rispetto a n, 
    % NB Va sommata perché ci sarebbe un meno moltiplicato a tutto che non ho messo
    dRp = diag(dt*M*dR_dp(n,p)); 

    J11 = A;
    J12 = M;
    J13 = -M;
    J21 = zeri;
    J22 = M - dRn;
    J23 = -dRp;
    J31 = zeri; 
    J32 = -dRn;
    J33 = M - dRp;
    jac = [J11, J12, J13;
           J21, J22, J23;
           J31, J32, J33];

end

end
