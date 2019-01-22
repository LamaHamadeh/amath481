% Forward Euler method

% Define tstart, tend and deltat
tstart = 0;
tend = 10;
dt = 0.3;

% set t vector
tspan = [tstart:dt:tend];

% predefine y vector
y = zeros(1,length(tspan));

% set initial condition
y(1) = 1;

% set the coefficient lambda of the linear ODE
lambda = -10;

% lambda = 10;
% lambda = -10;

% make the forward Euler iteration
for n=1:length(y)-1    
    y(n+1)=y(n) + dt*ffunc1(tspan(n),y(n),lambda);
end


%show the resulting y(t)
figure;

set(gca,'FontSize',16);
plot(tspan,y,'-o'); 
hold on;

% compare with the analitic solution if known
tspanfine=linspace(tstart,tend,10000); 
plot(tspanfine,exp(lambda*tspanfine),'LineWidth',3);  %exact solution x=exp(t)

xlabel('t','Fontsize',20); 
ylabel('y(t)','Fontsize',20); 
legend('Euler Approx','Exact Solution')