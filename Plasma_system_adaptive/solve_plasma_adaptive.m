options = optimoptions('fsolve','SpecifyObjectiveGradient',true); %'MaxIterations',800,'OptimalityTolerance',1e2,'StepTolerance',1e-6);

x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],1);
vbc_fun = @(tf) (Vendin - Vsrtin)/Tin * tf;


t = tsavein(1);
dt = dt0in;
breakflag = false;

% Loops over K
for it = 2 : numel (tsavein)  
    % Make sure that it executes exact times steps of tsave
    while t < tsavein(it)
        % clean leftovers
        if (t + dt > tsavein(it))
	        dt = tsavein(it) - t;
        end

        % + dt
        fun = @(x) Funz_plasma_adaptive(x,x_prec,[vbc_fun(t+dt);v_bcin(2)],n_bcin,p_bcin,dt,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = x_prec;
        [xnew, ~, ~, ~]  = fsolve(fun,initial_guess, options);
    
        % + dt/2
        funTemp1 = @(x) Funz_plasma_adaptive(x,x_prec,[vbc_fun(t+dt/2);v_bcin(2)],n_bcin,p_bcin,dt/2,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = x_prec;
        [xtemp, ~, exitflag, ~]  = fsolve(funTemp1,initial_guess,options);

        if exitflag == 0
            disp('Iteration limit reached, solution may not be accurate.');
            breakflag = true;
            break;
        elseif exitflag < 0
            disp('fsolve failed to find a solution.');
            breakflag = true;
            break;
        end

        funTemp2 = @(x) Funz_plasma_adaptive(x,xtemp,[vbc_fun(t+dt);v_bcin(2)],n_bcin,p_bcin,dt/2,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = xtemp;
        [xtemp, fval, exitflag, output] = fsolve(funTemp2,initial_guess, options);
        % optimset('TolX', tol/10, 'TolFun', tol/10))
        if (norm (xnew - xtemp) < tol)
            t = t+dt;
            dt = dt*1.2;
            fprintf("Difference smaller than tol\n")
            fprintf ('t = %g, dt = %g,  V(0) = %g\n', t*tbar, dt*tbar, vbc_fun(t)*Vbar)
            % pause(2)
            x_prec = xnew;
        else
            fprintf("FAIL, halving dt\n")
            fprintf("Norm = %g", norm (xnew - xtemp))
            % pause(2)
	        dt = dt/2;
        end
    end
    if breakflag
        break
    end

    fprintf("Saved solution n = %g, at time = %g, V(0) = %g\n", it, tsave(it), vbc_fun(t)*Vbar)
    pause(2);
    X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],it) = x_prec;
    X(1,it) = vbc_fun(t);

end

% for tol = 1e-3 dt hovers around ~7e-6 