function [lap] = laplacian(n)

    dx = 1/n;
    dy = 1/n;

    x = 0:dx:1;
    y = 0:dy:1;

    [X Y] = meshgrid(x,y) ;
    
    % Building the 1D Dirichelet minus Laplacian matrix
    e = ones(n-1,1);
    Asp = spdiags([-e 2*e -e], -1:1, n-1, n-1);
  
    % Building the 1D Neumann minus Laplacian matrix
    e = ones(n+2,1);
    Bsp = spdiags([-e 2*e -e], -1:1, n+1, n+1);
    Bsp(1,1:3) = [1.5 -2 0.5 ]*dy;
    Bsp(end,end-2:end) = [ 0.5 -2 1.5 ]*dy;
    
    % creating the identities (here carefull with the 
    % boundaries)
    I_A = speye(n+1,n+1);   % I_A(1,1) = 0;   I_A(end,end) = 0; not understanding the need for this
    I_B = speye(n-1,n-1);
    
    % assembling the 2D minus Laplacian
    lap = kron(Asp/dx^2,I_A) + kron(I_B,Bsp/dy^2);
    
    % writing the source term
    f = cos(2*pi*y);    % f(1) = 0;   f(end) = 0; likewise leaving this modificaiton out
    e_1 = zeros(n-1,1);   e_1(1) = 1;
    F = kron(e_1, f).';
    f = F(:)/dx^2;

end


