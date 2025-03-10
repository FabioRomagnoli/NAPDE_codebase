%% risolto con fsolve
for k=1:K

x_prec=X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);

% non-constant generation temr
fun = @(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, betain, Eiin);

% constant generation term
% fun = @(x) Funz_jacob_plasma_gen(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, gen);

initial_guess=x_prec; % lo prendo così,,,,,,,,)

[x,fval,exit,output,jacobian]=fsolve(fun,initial_guess,options);

% if k > 100
%     figure;
%     Ar_full = Plasma_rhs2(r, [v_bc(1,k); x(1:lr-2); v_bc(2,k)], munin, alphain, Vthin, -1);
%     Ar = Ar_full(2:end-1,2:end-1);
%     plot(r(2:end-1),Ar*x(lr-1:2*lr-4))
% 
%     figure;
%     plot(r(2:end-1),x(1:lr-2));
%     grid on;
% end

X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1)=x;

end

