solveType = " fsolve";
 

% risolto con fsolve
for k=1:K

x_prec=X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);

fun=@(x) Funz_Jacob(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,xin,Nin,muin,epsin,niin,Vthin,tauin);

initial_guess=x_prec; % lo prendo così,,,,,,,,)

% Nothing provided
% options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');

% Jacobian Pattern
% options = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e5,'MaxIterations',800,'JacobPattern',Jp,'OptimalityTolerance',1e2,'StepTolerance',1e-20);

% Full Jacobian
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);

[x,fval,exit,output,jacobian]=fsolve(fun,initial_guess,options);

X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1)=x;

end



















