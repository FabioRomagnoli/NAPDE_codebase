function F = ax_gradient(r)
% central_diff_matrix builds a sparse differentiation matrix F for a
% nonuniform 1D mesh r using central differences for interior points
% and one-sided differences for the endpoints.
%
% INPUT:
%   r - column vector of mesh points (assumed sorted in increasing order)
%
% OUTPUT:
%   F - sparse differentiation matrix such that F*v approximates the 
%       derivative of v with respect to r

    N = length(r);
    dr = diff(r);  % differences between consecutive mesh points, length N-1

    % For interior nodes i = 2:N-1, define h_im1 = r(i)-r(i-1) and h_i = r(i+1)-r(i)
    h_im1 = dr(1:end-1);  % for i = 2 to N-1, length = N-2
    h_i   = dr(2:end);    % for i = 2 to N-1, length = N-2

    % Compute the coefficients for the interior nodes based on quadratic interpolation:
    % f'(r_i) ~ c1*v(r_{i-1}) + c2*v(r_i) + c3*v(r_{i+1})
    L = - h_i ./ (h_im1 .* (h_i + h_im1));  % subdiagonal coefficients for rows 2:N-1
    M = (h_i - h_im1) ./ (h_i .* h_im1);      % main diagonal coefficients for rows 2:N-1
    U =   h_im1 ./ (h_i .* (h_i + h_im1));      % superdiagonal coefficients for rows 2:N-1

    % Build the interior part of F using spdiags.
    % This creates a sparse matrix of size (N-2) x N where:
    %  - the -1 diagonal (columns i-1) has entries L,
    %  - the 0 diagonal (columns i) has entries M,
    %  - the 1 diagonal (columns i+1) has entries U.
    F_interior = spdiags([L, M, U], [-1, 0, 1], N-2, N);

    % Allocate the full sparse differentiation matrix F
    F = spalloc(N, N, 3*N);
    F(2:N-1, :) = F_interior;

    % For the boundaries, use first-order differences:
    % Forward difference at the first node:
    h1 = dr(1);
    F(1,1) = -1/h1;
    F(1,2) =  1/h1;
    
    % Backward difference at the last node:
    hN = dr(end);
    F(N,N-1) = -1/hN;
    F(N,N)   =  1/hN;
end
