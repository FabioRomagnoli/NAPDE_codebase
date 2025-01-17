function [xvect,it] = newtonsys(x0, maxit, toll, fun, verbose) 

    dampit     = 10;
    dampcoeff  = 2;
        
    x_old = x0;
    
    [res,jac]=fun(x_old);
    
    normr   = norm (res, inf);
    normrnew   = normr;
    
    first_iteration_output (verbose, 0, normrnew);
    
    
    for it = 1:maxit
        delta = - jac \ res;
    
        tk = 1;         % Adaptive step
    
        for dit = 1:dampit
    
            x_new = x_old + tk * delta;
        
            [resnew,jacnew] = fun(x_new);
        
            normrnew = norm(resnew, inf);
        
            if (normrnew > normr)
                tk = tk / dampcoeff;
            else
                jac = jacnew;
                res = resnew;
                break
            end	
        end
    
        normx = norm(x_new - x_old,2);
        x_old = x_new;
        normr = normrnew;
        iteration_output (verbose, it, normrnew)
    
        if (normr <= toll || normx <= toll)
            break
        end
    end
    
    xvect = x_new;
end


function first_iteration_output (verbose, it, normrnew)
    if (verbose)
        fprintf ('Newton iteration  %3d', it);
        fprintf ('  Scaled system residual %6g\n', normrnew); 
    end
end

function iteration_output (verbose, it, normrnew)
    if (verbose)
        fprintf ('                  %3d', it);
        fprintf ('                         %6g\n', normrnew); 
    end
end
