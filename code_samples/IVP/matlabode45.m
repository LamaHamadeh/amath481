
% Define tstart, tend and deltat
tstart = 0;
tend = 10;
dt = 0.3;

% set t vector
tspan = [tstart:dt:tend];

% predefine y vector
y = zeros(1,length(tspan));

% set initial condition
y0 = 1;

% set the coefficient lambda of the linear ODE
lambda = 10;


%Call ODE solver ode45 - RK4

%[t,y]=ode45(@(t,y) ffunc1(t,y,lambda),[tstart,tend],y0);

% with tspan
%[t,y]=ode45(@(t,y) ffunc1(t,y,lambda),tspan,y0);

% with tolerance
% tol=10^-4;
% options = odeset('RelTol',tol,'AbsTol',tol);

%[t,y]=ode45(@(t,y) ffunc1(t,y,lambda),tspan,y0,options);

% with vector valued function

[t,y]=ode45(@(t,y) ffunc2(t,y,lambda),tspan,[1;1]);

y1 = y(:,1);
y2 = y(:,2);

%show the resulting y(t)
figure;

set(gca,'FontSize',16);
plot(t,y1,'-o'); 
hold on;

% compare with the analytic solution if known
tspanfine=linspace(tstart,tend,10000); 
plot(tspanfine,exp(lambda*tspanfine),'LineWidth',3);  %exact solution x=exp(t)

xlabel('t','Fontsize',20); 
ylabel('y(t)','Fontsize',20); 
legend('ODE45 Approx','Exact Solution')