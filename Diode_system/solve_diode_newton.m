% risolto con Newton
solveType = " Newton";

maxit = 150; 
toll = 1e-5;

for k = 1:K
   x_prec = X([2:lr-1 lr+2:2*lr-1 2*lr+2:end-1],k);

   fun=@(x) Funz_Jacob(x,x_prec,v_bcin(:,k),n_bcin(:,k),p_bcin(:,k),dtin,xin,Nin,muin,epsin,niin,Vthin,tauin);

   initial_guess = x_prec;  % E' ridotto

   [x_sol,it] = newtonsys(initial_guess, maxit, toll, fun, 1);

   X([2:lr-1 lr+2:2*lr-1  2*lr+2:3*lr-1],k+1) = x_sol;
end




















