%% AMATH 481 - PS 4 - Spencer Pease
%

clear all; close all;


%% Setup
%

%%% Define problem parameters %%%

tspan = 0:0.5:4;
nu = 0.001;
L=20;
n=64;


%%% Define grid space %%%

% problem space
x2 = linspace(-L/2, L/2, n+1); x = x2(1:n);
y=x;
[X, Y] = meshgrid(x, y);

% fourier space
kx = (2*pi/L)*[0:(n/2-1) (-n/2):-1]; kx(1)=10^-6;
ky=kx;
[KX, KY] = meshgrid(kx, ky);


%%% Setup basic RHS function %%%

[A, B, C] = createMatrices(n, L);
A(1,1) = 2 * (L/ n^2);

rhsFun = @(psi_vec, w_vec) ...
    (-B*psi_vec).*(C*w_vec) + (C*psi_vec).*(B*w_vec) + nu.*A*w_vec;


%%% Define initial conditions %%%

w0 = exp(-(X).^2 - ((Y.^2) / 20));
w0_vec = reshape(w0, n^2, 1);

%% Solve 
%

% options = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

%%% A\b %%%

tic;

[t1, w1Sol_vec] = ode45(...
    @(t, w_vec) vorticityBackslash(t, w_vec, rhsFun, A), ...
    tspan, w0_vec);

time1 = toc;
    

%%% LU %%%

[LL, U, P] = lu(A);

tic;

[t2, w2Sol_vec] = ode45(...
    @(t, w_vec) vorticityLU(t, w_vec, rhsFun, LL, U, P), ...
    tspan, w0_vec);

time2 = toc;


%%% FFT %%%

K = KX.^2 + KY.^2;

tic;

[t3, w3Sol_vec] = ode45(...
    @(t, w_vec) vorticityFFT(t, w_vec, rhsFun, K, n), ...
    tspan, w0_vec);

time3 = toc;


% Answers ---------------------------------------------------------------------

A1 = w1Sol_vec;
A2 = w2Sol_vec;
A3 = w3Sol_vec;


%% Plot and Fit
%

% wPlot_vec = w1Sol_vec;
% 
% for j = 1:size(wPlot_vec, 1)
%     
%     w_curr = reshape(wPlot_vec(j, :), n, n);
%     pcolor(X, Y, w_curr);
%     pause(0.01);
%     
% end


%% Write Data
%

save A1.dat  A1  -ascii
save A2.dat  A2  -ascii
save A3.dat  A3  -ascii


