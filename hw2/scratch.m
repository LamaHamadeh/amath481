
clear all; % clear all previously defined variables 
close all; % clear all previously defined figures


tol = 10^(-4); % define a tolerance level to be achieved
               % by the shooting algorithm 
col = ['r','b','g','c','m','k']; % eigenfunction colors

n0 = 100;     % define the parameter n0 
A = 1;        % define the initial slope at x=-1 
x0 = [0 A];   % initial conditions: x1(-1)=0, x1'(-1)=A 
xp = [-1 1];  % define the span of the computational domain


beta_start = n0; % beginning value of beta

for modes = 1:5 % begin mode loop
    
    beta = beta_start; % initial value of eigenvalue beta
    dbeta = n0/100; % default step size in beta
    
    for j = 1:1000 % begin convergence loop for beta
        
        [t,y] = ode45('shoot2', xp, x0, [], n0, beta); % solve ODEs
        
        if abs(y(end,1) - 0) < tol % check for convergence
            beta % write out eigenvalue
            break % get out of convergence loop
        end
        
        if (-1)^(modes+1)*y(end,1)>0 % this IF statement block
            beta = beta - dbeta;     % checks to see if beta
        else                         % needs to be higher or lower
            beta = beta + dbeta/2;   % and uses bisection to
            dbeta = dbeta/2;         % converge to the solution
        end                          %
        
    end % end convergence loop
    
    beta_start = beta - 0.1; % after finding eigenvalue, pick
                             % new starting value for next mode
    norm = trapz(t, y(:,1).*y(:,1)) % calculate the normalization
    plot(t, y(:,1)/sqrt(norm), col(modes)); hold on % plot modes
    
end % end mode loop
