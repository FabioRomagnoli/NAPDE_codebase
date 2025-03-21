function J = Comp_current(r,mu,q,v,Vth,z,n)
% X deve essere il vettore dei pt medi degli intervalli, o comunque deve
% essere un pt per ogni sottointerv per funzionare
  DV = z .* diff (v) / Vth;

  [Bp, Bn] = bimu_bernoulli (DV);

  nBn=n(2:end).*Bn;
  nBp=n(1:end-1).*Bp;

  J=2*pi*z*q*mu*Vth./(log(r(1:end-1) ./ r(2:end))).*(nBn-nBp);
  % units of J are A/m
end








