%% risolto con fsolve
for k=1:K

x_prec=X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);

fun = @(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, gen);

initial_guess=x_prec; % lo prendo così,,,,,,,,)

[x,fval,exit,output,jacobian]=fsolve(fun,initial_guess,options);

X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1)=x;

end

