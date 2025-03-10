function Ar = Plasma_rhs3(r, v, mu, alpha, Vth, z)
    ndof = numel(r);
    dr = diff(r);
    coeff = (dr/2) ./ log(r(1:end-1)./r(2:end));

    DV = z .* diff(v) / Vth;

    [Bp, Bm] = bimu_bernoulli(DV);

    fll = - (Bp .* coeff);   % fll corresponds to -B_{p_l} * coeff
    fru =   (Bm .* coeff);      % fru corresponds to B_{n_u} * coeff
    
    dd = [fll; 0] + [0; fru];
    du = [0; fru];    % upper diagonal (shifted down by one)
    dl = [fll; 0];   % lower diagonal (shifted up by one)
    
    Ar = - mu * Vth * alpha * spdiags([dl, dd, du], -1:1, ndof, ndof);

end