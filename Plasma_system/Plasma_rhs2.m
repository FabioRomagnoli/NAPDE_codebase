function Ar=Plasma_rhs2(r, v, mu, alpha, Vth, z)

lr = length(r);
dr=diff(r);

DV = z .* diff (v) / Vth;
[Bp, Bn] = bimu_bernoulli (DV);

%Scarta primo elem sopradiag e ultimo sottodiag

ld = [dr/2.*Bp./log(r(1:end-1)./r(2:end)); 0];
ud = [0; dr/2.*Bn./log(r(1:end-1)./r(2:end))];
dd = -ld-ud;
% Dovrebbe avere un modulo, per questo ho aggiunto un -
Ar = mu*Vth*alpha*spdiags([ld, dd, ud], -1:1, lr, lr);

end