function R=Plasma_rhs(r, v, n, mu, alpha, Vth, z)

dr=diff(r);


  DV = z .* diff (v) / Vth;
  [Bp, Bn] = bimu_bernoulli (DV);


nBp=n(1:end-1).*Bp;
nBn=n(2:end).*Bn;

xJ_l=abs(mu*Vth./(log(r(1:end-1) ./ r(2:end))).*(nBn-nBp).*dr/2);


xJ_u=abs(mu*Vth./(log(r(1:end-1) ./ r(2:end))).*(nBn-nBp).*dr/2);

R=alpha*([0; xJ_l]+[xJ_u; 0]);








