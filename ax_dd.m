function F = ax_dd (r, v, mu, Vth, z)

  ndof = numel (r);
  %h    = diff (r);
  
  DV = z .* diff (v) / Vth;
  [Bp, Bm] = bimu_bernoulli (DV);

  fll = - mu * Vth .* Bp ./ (log(r(1:end-1) ./ r(2:end)));
  fld =   mu * Vth .* Bm ./ (log(r(1:end-1) ./ r(2:end)));
  fru =   mu * Vth .* Bm ./ (log(r(1:end-1) ./ r(2:end)));
  frd = - mu * Vth .* Bp ./ (log(r(1:end-1) ./ r(2:end)));


  dd = [frd; 0] - [0; fld];
  du = [0; fru];
  dl = [-fll; 0];
  F = spdiags([dl, dd, du], -1:1, ndof, ndof);
  
end


%!demo
%! ua = .1; ub = 0; a = .7e-3; b = a + .15;
%! mu = .1; Vth = 26e-3; beta = 1;
%! msh = CreateTanhMshAsym (100, a, b-a, 1.001);
%! x = [a; msh.x(:); b];
%! phi =  (ua * log(b) - ub * log(a)) ./ log(b/a) + (ub - ua) * log(x) ./ log(b/a);
%! F = ax_dd (x, phi, .1, 26e-3, 1);
%! n = 1e18 * ones (size (x));
%! n(1) = 1e14;
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * n([1 end]));
%! semilogy (x, n, 'x-')

%!test
%! %% - d/dx (x mu Vth (n' - n phi' / Vth)) = x g(x)
%! a   = .7e-3; b = a + .15; 
%! mu  = .1; Vth = 26e-3; beta = 100;
%! msh = CreateTanhMsh(400, a, b-a, 2.4577e-5, -0.1);
%! x = [a; msh.x(:); b];
%! phi = - beta * x;
%! g   = -(Vth - beta) * mu * exp(-x) .* (-1 + x) ./ x;
%! F   = ax_dd (x, phi, mu, Vth, -1);
%! M   = ax_mass (x, 1);
%! b   = M*g;
%! n = nex = exp(-x);
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (b(2:end-1)-F(2:end-1, [1 end]) * n([1 end]));
%! assert (norm (n - nex, inf) < 1e-3)
