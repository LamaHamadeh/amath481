%% AMATH 481 - PS 1 - Spencer Pease
%

clear all; close all;


%% Exercise 1 - Euler and Heun
%

% Setup -----------------------------------------------------------------------

% Consider the ODE:
dydt = @(t, y) -3*y*sin(t);

% with the exact solution:
sol = @(t) pi*exp(3*(cos(t) - 1))/sqrt(2);

% and initial condition:
y0 = pi / sqrt(2);

% over a time interval:
t_int = [0, 5];

% with dela Ts:
dt = 2.^(-1*(2:8));

% Store error of methods for different time steps:
E_euler = zeros(1, length(dt));
E_heun = zeros(1, length(dt));


% Solve -----------------------------------------------------------------------

for j = 1:length(dt)
    
    y_true = sol(t_int(1):dt(j):t_int(2));
    y_num_euler = forwardEuler(dydt, y0, t_int, dt(j));
    y_num_heun = heun(dydt, y0, t_int, dt(j));
    
    E_euler(j) = mean(abs(y_true' - y_num_euler));
    E_heun(j) = mean(abs(y_true' - y_num_heun));
    
end


% Plot and Fit ----------------------------------------------------------------

% hold on;
% loglog(dt, E_euler, '-*b')
% loglog(dt, E_heun, '-or')
% title("Error vs dt")
% legend("Euler", "Heun")
% xlabel("Delta t")
% ylabel("Error")
% hold off;

fit_euler = polyfit(log(dt), log(E_euler), 1);
fit_heun  = polyfit(log(dt), log(E_heun), 1);


% Answers ---------------------------------------------------------------------

A1 = y_num_euler;
A2 = E_euler;
A3 = fit_euler(1);
A4 = y_num_heun;
A5 = E_heun;
A6 = fit_heun(1);


%% Exercise 2 - van der Pol
%

% Setup -----------------------------------------------------------------------

% Consider the van der Pol oscillator:
vdPol = @(t, y, eps) [y(2); -1*eps*(y(1)^2 - 1)*y(2) - y(1)];


% Solve -----------------------------------------------------------------------

% a) For the conditions
y0 = [sqrt(3); 1];
eps = [0.1, 1, 20];
t_span = 0:0.5:32;

Y = zeros(length(t_span), length(eps));

for j = 1:length(eps)
    
    [~, y_tmp] = ode45(@(t, y) vdPol(t, y, eps(j)), t_span, y0);
    Y(:, j) = y_tmp(:, 1);
    
end

% b) For the conditions
y0 = [2; pi^2];
eps = 1;
t_span = [0, 32];
TOL = 10.^(-1*(4:10));

method = {@ode45, @ode23, @ode113};
avg_dt = zeros(length(method), length(TOL));

for k = 1:length(method)
    
    for j = 1:length(TOL)
    
    options = odeset('AbsTol', TOL(j), 'RelTol', TOL(j));
    [T, ~] = method{k}(@(t, y) vdPol(t, y, eps), t_span, y0, options);
    avg_dt(k, j) = mean(diff(T));
    
    end
end


% Plot and Fit ----------------------------------------------------------------

fit_ode45  = polyfit(log(avg_dt(1, :)), log(TOL), 1);
fit_ode23  = polyfit(log(avg_dt(2, :)), log(TOL), 1);
fit_ode113 = polyfit(log(avg_dt(3, :)), log(TOL), 1);


% Answers ---------------------------------------------------------------------

A7  = Y;
A8  = fit_ode45(1);
A9  = fit_ode23(1);
A10 = fit_ode113(1);


%% Exercise 3 - Fitzhugh 
%

% Setup -----------------------------------------------------------------------

% Consider two Fitzhugh neurons:
a1 = 0.05; a2 = 0.25; b = 0.01; c = 0.01; I = 0.1;
fitzhugh = @(t, y, d12, d21) [
    -1*y(1)^3 + (1 + a1)*y(1)^2 - a1*y(1) - y(3) + I + d12*y(2);
    -1*y(2)^3 + (1 + a2)*y(2)^2 - a2*y(2) - y(4) + I + d21*y(1);
    b*y(1) - c*y(3);
    b*y(2) - c*y(4);
    ];

% with the initial conditions
y0 = [0.1, 0.1, 0, 0];
t_span = 0:0.5:100;


% Solve -----------------------------------------------------------------------

d12 = [0, 0, -0.1,-0.3, -0.5];
d21 = [0, 0.2, 0.2, 0.2, 0.2];

Y = zeros(length(t_span), length(y0), length(d12));

for j = 1:length(d12)
    
    [~, y_tmp] = ode15s(@(t,y) fitzhugh(t, y, d12(j), d21(j)), t_span, y0);
    Y(:,:,j) = y_tmp;
    
end


% Answers ---------------------------------------------------------------------

A11 = Y(:,:,1);
A12 = Y(:,:,2);
A13 = Y(:,:,3);
A14 = Y(:,:,4);
A15 = Y(:,:,5);


%% Write Data
%

save A1.dat  A1  -ascii
save A2.dat  A2  -ascii
save A3.dat  A3  -ascii
save A4.dat  A4  -ascii
save A5.dat  A5  -ascii
save A6.dat  A6  -ascii
save A7.dat  A7  -ascii
save A8.dat  A8  -ascii
save A9.dat  A9  -ascii
save A10.dat A10 -ascii
save A11.dat A11 -ascii
save A12.dat A12 -ascii
save A13.dat A13 -ascii
save A14.dat A14 -ascii
save A15.dat A15 -ascii
