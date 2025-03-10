figure;
title('Electron and hole concentrations')

for k=1:K
v = X(1:lr,k);
n = X(lr + 1:2*lr,k);
p = X(2*lr + 1:end,k);

% Descaling procedure
v = v*Vbar;
n = n*nbar;
p = p*nbar;


% semilogy(r,p,'o-',r,n,'x-')
semilogy(r,p,r,n)
legend("p","n")
grid on;
drawnow;
% pause(0.1)
end

% figure;
% plot(r,v)
% title('Electric potential')
% 
% 
% figure;
% phi00 = Vend * (1 - log (r/r0)/log(r1/r0));
% plot(r,v,r,phi00);
% title('Electric potential')
% 

figure;
x_medi = (r(1:end-1)+r(2:end))/2;
title('Corrente elettrica')
Jn = Comp_current(r,mun,q,v,Vth,-1,n);
Jp = Comp_current(r,mup,q,v,Vth, 1,p);
JJ = Jn+Jp;
hold on;

plot(x_medi,Jn,"-o");
plot(x_medi,Jp);
plot(x_medi,JJ);
legend('Jn','Jp','JJ')
hold off;