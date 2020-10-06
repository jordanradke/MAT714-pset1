%%% multigrid method 

% inputs: initial guess and source term
function [w] = multigrid(u, source)

n = length(u)-1;

%% presmoothing
% number of gauss-seidel iterations
nu = 3;

u_gs = gauss_seidel(nu, u, source);

%% compute residual
% source term
y = 0:(1/n):1;
f_bdy = cos(2*pi*y);    % f(1) = 0;   f(end) = 0; leaving this modificaiton out
e_1 = zeros(n-1,1);   e_1(1) = 1;
F = kron(e_1, f_bdy).';
f = n^2*F(:);
    
Delta = full(laplacian(n));  % laplacian

r = f - Delta*u_gs(:);


%% coarsen residual
r_half = restriction(r);


%% solve for coarse error. Direct solve if small enough
if n < 5
    Delta_half = full(laplacian(n-1));
    e_half = Delta_half\r_half;
else
    % main recursive step
    r_half = reshape(r_half, [n,n-2]);
    e_half_guess = zeros(n,n-2);
    multigrid(e_half_guess,-r_half);
end

%% interpolate e back to original grid, use it to update u_gs
e_inter = interpolation(e_half, (n+1)*(n-1));
e = reshape(e_inter, [n+1,n-1]);
u_new = u_gs - e;


%% post-smoothing: use u_new as your starting guess for GS
w = gauss_seidel(nu, u_new, homog_source);

end





