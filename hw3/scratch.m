
m = 8;                  % N value in x and y
n = m^2;                % total matrix size

e0 = zeros(n, 1);       % base vector of zeros
e1 = ones(n, 1);        % base vector of ones


e22 = e1; % copy the one vector 
e44 = e0; % copy the zero vector 

for j = 1:m 
    e22(m*j) = 0; % overwrite every m^th value with zero 
    e44(m*j) = 1; % overwirte every m^th value with one 
end

e33(2:n, 1) = e22(1:n-1, 1); e33(1, 1) = e22(n, 1); % shift to correct
e55(2:n, 1) = e44(1:n-1, 1); e55(1, 1) = e44(n, 1); % positions


matA = spdiags([e1 e1 e55 e22 -4*e1 e33 e44 e1 e1], [-(n-m) -m -m+1 -1 0 1 m-1 m (n-m)],n,n);