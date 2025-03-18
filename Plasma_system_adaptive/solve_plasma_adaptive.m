x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],1);

vbc_fun = @(dtv) (Vendin - Vsrtin)/Tin * dtv + v_bcin(1);


t = tsavein(1);
dt = dt0in;

% Loops over K
for it = 2 : numel (tsavein)  
    % Make sure that it executes exact times steps of tsave
    while t < tsavein(it)
        % clean leftovers
        if (t + dt > tsavein(it))
	        dt = tsavein(it) - t;
        end

        % + dt
        fun = @(x) Funz_plasma_adaptive(x,x_prec,[vbc_fun(dt);v_bcin(2)],n_bcin,p_bcin,dt,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = x_prec;
        [xnew, ~, ~, ~]  = fsolve(fun,initial_guess);
    
        % + dt/2
        funTemp1 = @(x) Funz_plasma_adaptive(x,x_prec,[vbc_fun(dt/2);v_bcin(2)],n_bcin,p_bcin,dt/2,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = x_prec;
        [xtemp, ~, exitflag, ~]  = fsolve(funTemp1,initial_guess);


        if exitflag == 0
            disp('Iteration limit reached, solution may not be accurate.');
            break;
        elseif exitflag < 0
            disp('fsolve failed to find a solution.');
            break;
        end


        funTemp2 = @(x) Funz_plasma_adaptive(x,xtemp,[vbc_fun(dt);v_bcin(2)],n_bcin,p_bcin,dt/2,rin,munin,mupin,epsin,alphain,Sin,Vthin);
        initial_guess = xtemp;
        [xtemp, fval, exitflag, output] = fsolve(funTemp2,initial_guess);
        % ,optimset('TolX', tol/10, 'TolFun', tol/10)
 
        if (norm (xnew - xtemp) < tol)
            v_bcin(1) = vbc_fun(dt);
            t = t+dt;
            dt = dt*1.2;
            fprintf("Difference smaller then tol\n")
            fprintf ('t = %g, dt = %g\n', t*tbar, dt*tbar)
            % pause(2)
            x_prec = xnew;
        else
            fprintf("FAIL, halfing dt\n")
            fprintf("Norm = %g", norm (xnew - xtemp))
            % pause(2)
	        dt = dt/2;
        end
    end
    fprintf("Saved solution n = %g at time = %g", it, tsave(it))
    pause(2);
    X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],it) = x_prec;
    X(1,it) = v_bcin(1);

end

