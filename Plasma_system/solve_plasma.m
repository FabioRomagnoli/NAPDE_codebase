%% risolto con fsolve
for k=1:K    
    x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);
    
    if solveConstGen
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','iter'); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);
        fun = @(x) Funz_jacob_plasma_gen(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, genin);
    else 
        options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','Display','iter');
        fun = @(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin,betain,Eiin);
    end


    initial_guess = x_prec; % lo prendo così,,,,,,,,)


    [x,fval,exit,output,jacobian] = fsolve(fun,initial_guess,options);
    
    X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1) = x;
    
end

