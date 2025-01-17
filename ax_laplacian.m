function F = ax_laplacian (r, mu)

  ndof = numel (r);
  fl = mu ./ log(r(1:end-1) ./ r(2:end));
  dd = [-fl; 0] + [0; -fl];
  du = [0; fl];
  dl = [fl; 0];
  F = spdiags([dl, dd, du], -1:1, ndof, ndof);
 
end

%!demo
%! ua = 100; ub = 0; a = .7e-3; b = a + .15;
%! msh = CreateTanhMshAsym (100, a, b-a, 1.001);
%! x = [a; msh.x(:); b];
%! F = ax_laplacian (x, 1);
%! n = zeros (size (x));
%! n(1) = 100;
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * n([1 end]));
%! phi =  (ua * log(b) - ub * log(a)) ./ log(b/a) + (ub - ua) * log(x) ./ log(b/a);
%! plot (x, n, 'bx', x, phi)

%!test
%! a = .7e-3; b = a + .15;
%! mu = .1; Vth = 26e-3; 
%! msh = CreateTanhMshAsym (100, a, b-a, 1.001);
%! x = [a; msh.x(:); b]; 
%! g = -(mu*Vth*exp(-x).*(-1 + x)) ./ x;
%! F = ax_laplacian (x, mu*Vth);
%! M = ax_mass (x, 1);
%! b = M*g;
%! n = nex = exp(-x);
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (b(2:end-1)-F(2:end-1, [1 end]) * n([1 end]));
%! assert (norm (n - nex, inf) <= 1e-4)

%!test
%! a = .7e-3; b = a + .15; 
%! mu = .1; Vth = 26e-3; 
%! msh = CreateTanhMshAsym (100, a, b-a, 1.001);
%! x = [a; msh.x(:); b]; 
%! F = ax_laplacian (x, mu*Vth);
%! n = nex = log(x/b) ./ log(a/b);
%! n(2:end-1) = F(2:end-1, 2:end-1) \ (-F(2:end-1, [1 end]) * n([1 end]));
%! assert (norm (n - nex, inf) <= 1e-10)
