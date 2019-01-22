function [A, B, C] = createMatrices(m, L)

% Setup -----------------------------------------------------------------------

% m = 8;                  % N value in x and y
n = m^2;                % total matrix size


% d = (10 - (-10)) / m;   % delta x, delta y
d = L / m;

% aScale = -4             % factor to scale the main diagonal of A

% e0 = zeros(n, 1);       % base vector of zeros
e1 = ones(n, 1);        % base vector of ones


% Solve -----------------------------------------------------------------------

%%% Create matrix A %%%

e2 = ones(n, 1); e2(mod(0:n-1, m) == 0) = 0; 
e3 = zeros(n, 1); e3(mod(1:n, m) == 0) = 1;

ind_a = [1 m-1 m (n-m)];
diag_a = [e2 e3 e1 e1];

A = (1 / d^2) * ...
    spdiags([rot90(diag_a, 2) -4*e1 diag_a], [-fliplr(ind_a) 0 ind_a], n, n);


%%% Create matrix B %%%

ind_b = [m n-m];
diag_b = [e1 -e1];

B = (1/(2*d)) * ...
    spdiags([-1*rot90(diag_b, 2), diag_b], [-fliplr(ind_b), ind_b], n, n);


%%% Create matrix C %%%

e4 = repmat([0; ones(m-1, 1)], m, 1);
e5 = repmat([zeros(m-1, 1); -1], m, 1);

ind_c = [1, m-1];
diag_c = [e4, e5];

C = (1/(2*d)) * ...
    spdiags([-1*rot90(diag_c, 2), diag_c], [-fliplr(ind_c), ind_c], n, n);

end
