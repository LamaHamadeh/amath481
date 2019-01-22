%%% main script %%%

close all;clc

tol = 10^-4;
n0=100;
xp = [-1 1];
A=1;
x0 = [0;A];
beta_start = n0;



for modes = 1:5
    
beta= beta_start;
dbeta = n0/100;


for j=1:1000
    
    [t,y] = ode45( @(t,y) rhsfunc(t,y,n0,beta),xp,x0);
    
    if abs(y(end,1)-0)<tol
        beta
        break;
    end
    
    if (-1)^(modes+1)*(y(end,1))>0
        beta = beta-dbeta;
        
    else
        beta = beta + dbeta/2;
        dbeta = dbeta/2;
    end
    
    
end

beta_start = beta-0.1;

plot(t,y(:,1)); hold on;

    
end



    
    
    
    
    
    
    
    
    
    