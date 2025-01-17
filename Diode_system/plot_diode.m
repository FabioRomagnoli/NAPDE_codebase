v = X(1:lr,end);
n = X(lr+1:2*lr,end);
p = X(2*lr+1:end,end);

% Descaling procedure
r = r*xbar;
v = v*Vbar;
n = n*nbar;
p = p*nbar;

figure
plot(r,v)
title("Electric potential" + solveType)

figure
semilogy(r,p,r,n,r, sqrt(n.*p))
title("Electron and hole concentrations" + solveType)
legend ('p', 'n','sqrt(n*p)')

