function F = ax_mass (r, mu)
  dd = [0; diff(r)/2] + [diff(r)/2; 0];
  dd = mu .* dd .* r;
  F = diag (dd);
end

%!test
%! a = .7e-3; b = a + .15;
%! mu = .1; Vth = 26e-3; 
%! msh = CreateTanhMshAsym (100, a, b-a, 1.001);
%! x = msh.x(:);
%! g = -(mu*Vth*exp(-x).*(-1 + x)) ./ r;
%! F = ax_laplacian (x, mu*Vth);
%! M = ax_mass (x, 1);
%! b = M*g;
%! n = nex = exp(-x);
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (b(2:end-1)-F(2:end-1, [1 end]) * n([1 end]));
%! assert (norm (n - nex, inf) <= 1e-4)
