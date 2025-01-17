solveType = " Split";


% Nothing provided
% options_reaction = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');
 
% Full Jacobian
options_reaction = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);

% Nothing provided
% options_transport = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');

% Full Jacobian
options_transport = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);


x_part=X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],1);

for k=1:K


% Reaction step
disp("reaction")
fun_reaction = @(x) reaction_fully_explicit(x,x_part,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,xin,Nin,muin,epsin,niin,Vthin,tauin);
[x_part,fval,exit,output,jacobian]=fsolve(fun_reaction,x_part,options_reaction);


% Transport step
disp("transport")
fun_transport=@(x) Funz_split_transport(x,x_part,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,xin,Nin,muin,epsin,niin,Vthin,tauin);
[x_part,fval,exit,output,jacobian]=fsolve(fun_transport,x_part,options_transport);


end

X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],end)=x_part;





















