% gauss-seidel with successive over-relaxation

k = 6;
n = 2^k;   % dyadic so points overlap 
h = 1/n;   % grid coarseness

bdy     = 0:h:1;
f       = sin(pi*bdy);       % bdy data @ x = 0
u       = zeros(n+1,n+1);     % initial guess
u_gs    = zeros(n+1,n+1);     % store gauss-seidel guesses here
maxiter = 100;              % what is reasonable here?
omega   = 2/(1+sin(pi*h));   % relaxation paramater

% remember: matlab indices start at 1 rather than 0.
for iter = 0:maxiter
    for i = 1:n+1                        % y variable
        for j = 1:n+1                    % x variable
            
            %% for the first row:
            if i == 1
               if j == 1
                   % dirichlet here
                   u(i,j) = f(i);
               elseif j < n
                   % neumann here
                   u_gs(i,j) = 4/3*u(i+1,j) - 1/3*u(i+2,j);
                   u(i,j) = u(i,j) + omega*(u_gs(i,j) - u(i,j));
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
                    u_gs(i,j) = .25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1));
                    u(i,j) = u(i,j) + omega*(u_gs(i,j) - u(i,j));

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
               elseif j < n
                   % neumann here
                   u_gs(i,j) = 4/3*u(i-1,j) - 1/3*u(i-2,j);
                   u(i,j) = u(i,j) + omega*(u_gs(i,j) - u(i,j));
               else
                   % dirichlet = 0 here.
                   u(i,j) = 0;  
               end
            end
            
        end
    end
end

% plot graph of u
surf(bdy,bdy,u)
   

%% comparison with 
        
            
        
            
