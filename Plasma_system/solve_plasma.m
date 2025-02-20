%% risolto con fsolve
for k=1:K

% Nothing provided
% options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');

% Full Jacobian
%options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);


x_prec=X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);

fun=@(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin);

initial_guess=x_prec; % lo prendo così,,,,,,,,)

[x,fval,exit,output]=fsolve(fun,initial_guess);

X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1)=x;

end

