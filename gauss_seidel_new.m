%%% gauss-seidel algorithm
function [w] = gauss_seidel_new(nu,u0, source)

n = length(u0)-1;
h = 1/n;   % grid coarseness

y       = 0:h:1;
f       = sign(cos(2*pi*y)).';         % bdy data @ x = 0
maxiter = nu;                  % what is reasonable here?
omega   = 2/(1+sin(pi*h));     % relaxation paramater

source = [zeros(n+1,1), source, zeros(n+1,1)];
u = [f, u0, zeros(n+1,1)];

% remember: matlab indices start at 1 rather than 0.
for iter = 0:maxiter
    for i = 1:n+1                        % y variable
        for j = 1:n+1                    % x variable
            
            %% for the first row:
            if i == 1
               if j == 1
                   % dirichlet here
                   u(i,j) = f(i);
               elseif j < n+1
                   % neumann here
                   u(i,j) = 4/3*u(i+1,j) - 1/3*u(i+2,j);
               else
                   % dirichlet = 0 here.
                   u(i,j) = 0;  
               end
            end
               
               
            %% for a typical row:
            if i > 1 && i < n+1
                % first entry in row is just left side boundary condition:
                if j == 1
                    u(i,j) = f(i);
                elseif j < n+1
                    u(i,j) = .25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1)) - .25*h^2*source(i,j);
                else
                    u(i,j) = 0;   % this should remain true since initialize zero
                                  % but just to be safe.
                end
            end
            
            %% for the top row:
            if i == n+1
               if j == 1
                   % dirichlet here
                   u(i,j) = f(i);
               elseif j < n+1
                   % neumann here
                   u(i,j) = 4/3*u(i-1,j) - 1/3*u(i-2,j);
               else
                   % dirichlet = 0 here.
                   u(i,j) = 0;  
               end
            end
            
        end
    end
end

w = u(:,2:n);

end





