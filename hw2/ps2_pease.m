%% AMATH 481 - PS 2 - Spencer Pease
%

clear all; close all;


%% Exercise 1 - Linear Shooting
%

% Setup -----------------------------------------------------------------------

TOL = 1e-4;
K = 1;
L = 4;
A = 1;              % phi'(-L) = A
bc = [1, 0];        % boundary conditions phi(-L) and phi(L)
y0 = [bc(1), A];    % initial conditions [phi(-L), phi(L)']
x_span = (-1*L):0.1:L;

max_mode = 5;
max_iter = 1000;

func = @(x, y, eps, K) [y(2); (K*x^2 - eps)*y(1)];


% Solve -----------------------------------------------------------------------

eigen_val = zeros(max_mode, 1);
eigen_fun = zeros(length(x_span), max_mode);

iters = ones(max_mode, 1).*-1;

eps_start = 8; %(K*L^2) / 2;

for modes = 1:max_mode
    
    eps = eps_start;
    deps = eps_start/100;
    
    for j = 1:max_iter 
        
        [t, y] = ode45(@(x, y) func(x, y, eps, K), x_span, y0);
        
        opt = y(end, 2) + sqrt(K*L^2 - eps)*y(end, 1); % value to optomize
        
        if abs(opt) < TOL
            iters(modes) = j;
            break
        end
        
        if (-1)^(modes+1)*opt > 0    % this IF statement block
            eps = eps - deps;        % checks to see if eps
        else                         % needs to be higher or lower
            eps = eps + deps/2;      % and uses bisection to
            deps = deps/2;           % converge to the solution
        end                          %
    end
    
    eigen_val(modes) = eps;
    eigen_fun(:, modes) = abs(y(:, 1) / sqrt(trapz(t, y(:, 1).^2)));
    
    eps_start = eps - .1;
    
end


% Plot and Fit ----------------------------------------------------------------

plot(t, eigen_fun);


% Answers ---------------------------------------------------------------------



%% Exercise 2 - van der Pol
%


%% Write Data
%


