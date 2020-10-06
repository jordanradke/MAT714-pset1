function [s] = restriction(r)

% create restriction matrix. The averaging at the boundaries is a little
% messy but hopefully this doesn't matter
n = length(r); 

% restriction matrix if even
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
   
s = R*r;

end