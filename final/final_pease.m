%% AMATH 481 - Final - Spencer Pease
%

clear all; close all;


%% Setup
%

%%% Define problem parameters

n = 16;
L = 2*pi;
tspan = 0:0.5:4;

A = [-1 -1 -1];
B = -A;

%%% Define problem space %%%

xyz = linspace(-L/2, L/2, n+1);
[X, Y, Z] = meshgrid(xyz(1:n), xyz(1:n), xyz(1:n));

x_vec = reshape(X, n^3, 1);
y_vec = reshape(Y, n^3, 1);
z_vec = reshape(Z, n^3, 1);


%%% Define fourier space %%%

kxyz = (2*pi/L)*[0:(n/2-1) (-n/2):-1]; kxyz(1) = 10^-6; 
[KX, KY, KZ] = meshgrid(kxyz, kxyz, kxyz);


%%% Define Laplacian %%%

K = KX.^2 + KY.^2 + KZ.^2;
Lap = -K;


%% Part A
%

% Setup -----------------------------------------------------------------------

psi0 = cos(X).*cos(Y).*cos(Z);
psi0f = fftn(psi0);
psi0f_vec = reshape(psi0f, n^3, 1);


% Solve -----------------------------------------------------------------------

[t1, psifSol1_vec] = ode45(...
    @(t, psif_vec) rhs(t, psif_vec, Lap, A, B, x_vec, y_vec, z_vec, n), ...
    tspan, psi0f_vec);


% Answers ---------------------------------------------------------------------

A1 = real(psifSol1_vec);
A2 = imag(psifSol1_vec);


%% Part B
%

% Setup -----------------------------------------------------------------------

psi0 = sin(X).*sin(Y).*sin(Z);
psi0f = fftn(psi0);
psi0f_vec = reshape(psi0f, n^3, 1);


% Solve -----------------------------------------------------------------------

[t2, psifSol2_vec] = ode45(...
    @(t, psif_vec) rhs(t, psif_vec, Lap, A, B, x_vec, y_vec, z_vec, n), ...
    tspan, psi0f_vec);


% Answers ---------------------------------------------------------------------

A3 = real(psifSol2_vec);
A4 = imag(psifSol2_vec);


%% Plot
%

psifPlot_vec = psifSol1_vec;

for j = 1:size(psifPlot_vec, 1)
    
    psifCurr_vec = psifPlot_vec(j, :);
    psifCurr = reshape(psifCurr_vec, n, n, n);
    psiCurr = ifftn(psifCurr);
    psiAbsCurr = psiCurr.*conj(psiCurr);
    
%     max(psiAbsCurr(:))
%     min(psiAbsCurr(:))

    isosurface(X, Y, Z, psiAbsCurr, .5)
    axis('square')
    pause(.05)
    
end


%% Write Data
%

save A1.dat  A1  -ascii
save A2.dat  A2  -ascii
save A3.dat  A3  -ascii
save A4.dat  A4  -ascii
