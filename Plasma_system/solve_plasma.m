%% risolto con fsolve
for k=1:K    
    x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);
    

    if false
        M_full = ax_mass(rin, 1);
        M = M_full(2:end-1,2:end-1);

        R = zeros(lr-2,1);
        R_full = Plasma_rhs(rin, [v_bcin(1); x_prec(1:lr-2); v_bcin(2)], [n_bcin(1); x_prec(lr-1:2*lr-4); n_bcin(2)], munin, alphain, Vthin, -1);
        % R = R_full(2:end-1);
        R(1:35) = R_full(2:36);

        % Ar_full = Plasma_rhs2(rin, [v_bcin(1,k); x_prec(1:lr-2); v_bcin(2,k)], munin, alphain, Vthin, -1);
        % Ar = Ar_full(2:end-1,2:end-1);
        % Ar_bc = Ar_full(2:end-1,[1 end]);
        % Rhs_red = zeros(lr-2,1);
        % nin = x_prec(lr-1:2*lr-4);
        % Rhs_red(1:36) = Ar(1:36,1:36)*nin(1:36);
        % Rhs_red(36) = Rhs_red(36)+Ar(36,37)*nin(37);

        % Rnew = M*genin + (R - M*genin)./2;
        
        figure();
        title('Generation Terms (adimensional) at previous step')

       
        hold on; 
        semilogy(rin(2:end-1), M*genin, "b-o", 'DisplayName', 'Const generation');
        semilogy(rin(2:end-1), R, "r-*", 'DisplayName', '\alpha*Jn');
        % semilogy(rin(2:end-1), Rnew, "k-s", 'DisplayName', 'rhs1New');
        % semilogy(rin(2:end-1), Rhs_red+Ar_bc(:,1)*n_bcin(1), "r-*", 'DisplayName', 'Rh2')
        hold off;

        legend();
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        grid on;
         
        % summDiff = sum(M*genin - R);
        % fprintf('summDiff = %.3f\n', summDiff);

        keyboard 
    end


    if solveConstGen
        options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','SpecifyObjectiveGradient',true,'Display','iter'); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);
        fun = @(x) Funz_jacob_plasma_gen(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin, genin);
    else 
        options = optimoptions('fsolve','Algorithm', 'trust-region-dogleg','Display','iter');
        % options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg','SpecifyObjectiveGradient',true,'Display','iter'); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);
        fun = @(x) Funz_jacob_plasma(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,rin,munin,mupin,epsin,alphain,Sin,Vthin,betain,Eiin, genin);
    end


    initial_guess = x_prec; % lo prendo così,,,,,,,,)

    [x,fval,exit,output,jacobian] = fsolve(fun,initial_guess,options);
    

    if false
        M_full = ax_mass(rin, 1);
        M = M_full(2:end-1,2:end-1);

        RSol = zeros(lr-2,1);
        R_full_sol = Plasma_rhs(rin, [v_bcin(1); x(1:lr-2); v_bcin(2)], [n_bcin(1); x(lr-1:2*lr-4); n_bcin(2)], munin, alphain, Vthin, -1);
        % R = R_full(2:end-1);
        RSol(1:35) = R_full_sol(2:36);
        
        % Ar_full = Plasma_rhs2(rin, [v_bcin(1,k); x(1:lr-2); v_bcin(2,k)], munin, alphain, Vthin, -1);
        % Ar = Ar_full(2:end-1,2:end-1);
        % Ar_bc = Ar_full(2:end-1,[1 end]);
        % Rhs_red = zeros(lr-2,1);
        % nin = x(lr-1:2*lr-4);
        % Rhs_red(1:36) = Ar(1:36,1:36)*nin(1:36);
        % Rhs_red(36) = Rhs_red(36)+Ar(36,37)*nin(37);

        % Rnew = M*genin + (RSol - M*genin)./2;
        figure();
        title('Generation Terms (adimensional) at subsequent step')
        hold on; 
        semilogy(rin(2:end-1), M*genin, "b-o", 'DisplayName', 'Const Generation');
        semilogy(rin(2:end-1), RSol, "r-*", 'DisplayName', '\alpha*Jn');

        % semilogy(rin(2:end-1), Rnew, "k-s", 'DisplayName', 'rhs1new');
        % semilogy(rin(2:end-1), Rhs_red+Ar_bc(:,1)*n_bcin(1), "r-*", 'DisplayName', 'Rh2')
        hold off; 

        legend();
        set(gca, 'YScale', 'log')
        % set(gca, 'XScale', 'log')

        grid on;
        keyboard  
    end


    if  false
        figure();
        vsol = x(1:lr-2);
        nsol = x(lr-1:2*lr-4);
        psol = x(2*lr-3:end);


        vprec = x_prec(1:lr-2);
        nprec = x_prec(lr-1:2*lr-4);
        pprec = x_prec(2*lr-3:end);


        hold on;
        semilogy(rin(2:end-1), nsol, "k-o", "DisplayName", "nsol");
        % semilogy(rin(2:end-1), psol, "k-o", "DisplayName", "psol" );
        semilogy(rin(2:end-1), nprec,"b-x", "DisplayName", "nprec");
        % semilogy(rin(2:end-1), pprec,"b-x", "DisplayName", "pprec");
        hold off;

        legend();
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

        grid on;
        keyboard
    end


    X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k+1) = x;
    
end

