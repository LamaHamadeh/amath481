%% AMATH 481 - PS 5 - Spencer Pease
%

clear all; close all;


%% Setup
%

L = 20;
m = 1; % number of spirals 
beta = 1;
D1 = 0.1;
D2 = 0.1;
tspan = 0:0.5:4;


%% FFT
%

% Setup -----------------------------------------------------------------------

n = 64;

xy = linspace(-L/2, L/2, n+1);
[X, Y] = meshgrid(xy(1:n), xy(1:n));

kxy = (2*pi/L)*[0:(n/2-1) (-n/2):-1]; kxy(1) = 10^-6; 
[KX, KY] = meshgrid(kxy, kxy);
K = KX.^2 + KY.^2;

fftLap = -K;

uvf0_vec = initialize_uvf(X, Y, m, n);


% Solve -----------------------------------------------------------------------

[t1, uvfSol_vec] =ode45(...
    @(t, uvf_vec) fft_rhs(t, uvf_vec, n, fftLap, D1, D2, beta), ...
    tspan, uvf0_vec);


% Answers ---------------------------------------------------------------------

A1 = real(uvfSol_vec);
A2 = imag(uvfSol_vec);


%% Chebyshev 
%

% Setup -----------------------------------------------------------------------

n = 31;

[D_unscaled, x_unscaled] = cheb(n-1);

D = (2/L).*D_unscaled;
x = (L/2).*x_unscaled;

DD = D^2;
DD(1, :) = zeros(1,n); DD(n, :) = zeros(1,n);

y = x;
[X, Y] = meshgrid(x, y);

uv0_vec = initialize_uv(X, Y, m, n);

I = eye(length(DD));
ChebLap = kron(DD, I) + kron(I, DD);


% Solve -----------------------------------------------------------------------

[t2, uvSolCheb_vec] =ode45(...
    @(t, uv_vec) cheb_rhs(t, uv_vec, n, ChebLap, D1, D2, beta), ...
    tspan, uv0_vec);


% Answers ---------------------------------------------------------------------

A3 = uvSolCheb_vec;


%% Testing

% uvPlot_vec = uvSolCheb_vec;
% % uvPlot_vec = uvfSol_vec;
% 
% for j = 1:size(uvPlot_vec, 1)
%     
%     uv_curr = uvPlot_vec(j, :);
%     
%     U_curr = reshape(uv_curr(1:n^2), n, n);
%     V_curr = reshape(uv_curr((n^2 + 1):end), n, n);
%     
% %     U_curr = real(ifft2(U_curr));
% %     V_curr = real(ifft2(V_curr));
%     
%     set(gcf, 'Position',  [600, 300, 800, 400]);
%     
%     ax1 = subplot(1, 2, 1);
%     pcolor(X, Y, U_curr);
%     title('U current')
%     axis('square')
%     
%     ax2 = subplot(1, 2, 2);
%     pcolor(X, Y, V_curr);
%     title('V current')
%     axis('square')
%     
%     pause(0.15);
%     
% end


%% Write Data
%

save A1.dat  A1  -ascii
save A2.dat  A2  -ascii
save A3.dat  A3  -ascii
