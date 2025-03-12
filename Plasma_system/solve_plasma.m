%% risolto con fsolve
for k=1:K

if k == 1
    x_prec = Xsol([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1]);
else
    x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);
end

%% Solver options:
% Nothing provided
options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');
% Full Jacobian
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter'); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);

%% Generation term: 
% non-constant generation term
fun = @(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin,betain,Eiin);
% constant generation term
% fun = @(x) Funz_jacob_plasma_gen(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, genin);

initial_guess=x_prec; % lo prendo così,,,,,,,,)


[x,fval,exit,output,jacobian]=fsolve(fun,initial_guess,options);


X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1)=x;

end

