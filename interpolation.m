function [s] = interpolation(r, n_new)

%% a messy loop to create restriction matrix
k = length(r);
n = n_new;   

% restriction matrix if n even
if mod(n,2) == 0
    
    R = zeros(n/2,n);
    
    for i = 1:n/2
        for j = 1:n
            % first row
            if i == 1
                if j == 2
                    R(i,j) = .5;
                    R(i,j-1) = .5;
                else
                    continue;
                end
            % typical row
            else
                if j == 2*i-1
                    R(i,j-1) = .25;
                    R(i,j)   = .5;
                    R(i,j+1) = .25;
                end
            end
        end
    end

% restriction matrix if n odd
else
    R = zeros((n+1)/2,n);
    
    for i = 1:(n+1)/2
        for j = 1:n
            % first row
            if i == 1
                if j == 2
                    R(i,j)   = .5;
                    R(i,j-1) = .5;
                else
                    continue;
                end
            % typical row
            elseif i < (n+1)/2
                if j == 2*i-1
                    R(i,j-1) = .25;
                    R(i,j)   = .5;
                    R(i,j+1) = .25;
                end
            % last row
            else
                if j == n
                    R(i,j)   = .5;
                    R(i,j-1) = .5;
                end
            end
        end
    end
end

%% interpolate r 
I = 2*R.';
I(2,1) = .5;
I(n,k) = 1;

s = I*r;

end