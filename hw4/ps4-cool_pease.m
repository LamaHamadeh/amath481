%% AMATH 481 - PS 4 - Spencer Pease
% As cool as it gets
%

clear all; close all;


%% Setup
%

%%% Define problem parameters %%%

tspan = 0:0.5:246;
nu = 0.001;
L=20;
n=128;


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


%% Initial Conditions
%

br = 5;
tl = -5;

w0 = ...
    exp(-1*(X+2.8 - br).^2 - 2*(Y + br).^2) - ...
    exp(-1*(X-2.8 - br).^2 - 2*(Y + br).^2) + ...
    exp(-1*(X - br).^2 - 2*(Y+2.8 + br).^2) - ...
    exp(-1*(X - br).^2 - 2*(Y-2.8 + br).^2) ...
    + 0 - ...
    exp(-1*(X+2.8 - tl).^2 - 2*(Y + tl).^2) + ...
    exp(-1*(X-2.8 - tl).^2 - 2*(Y + tl).^2) - ...
    exp(-1*(X - tl).^2 - 2*(Y+2.8 + tl).^2) + ...
    exp(-1*(X - tl).^2 - 2*(Y-2.8 + tl).^2);

w0 = w0 + rot90(-w0);

w0_vec = reshape(w0, n^2, 1);


%% Solve
%

K = KX.^2 + KY.^2;

[t4, wCool_vec] = ode45(...
    @(t, w_vec) vorticityFFT(t, w_vec, rhsFun, K, n), ...
    tspan, w0_vec);


%% Plot
%

for j = 1:size(wCool_vec, 1)
    
    w_curr = reshape(wCool_vec(j, :), n, n);
    
    set(gcf, 'Position',  [600, 100, 1201, 501]);
    
    ax1 = subplot(1, 2, 1);
    pcolor(X, Y, w_curr); 
    shading interp; colormap(lines(45)); axis off;
    
    ax2 = subplot(1, 2, 2);
    surf(X, Y, w_curr); 
    shading interp; colormap(lines(45)); axis off;
    zlim([-1, 1])
    
    
    pause(0.005);
    
    
end


%% Save images
%

mkdir('images');

for j = 1:size(wCool_vec, 1)
    
    w_curr = reshape(wCool_vec(j, :), n, n);
    
    figure(j);
    set(j, 'visible', 'off');
    
    set(gcf, 'Position',  [600, 100, 1201, 501]);
    
    ax1 = subplot(1, 2, 1);
    pcolor(X, Y, w_curr); 
    shading interp; colormap(lines(45)); axis off;
    
    ax2 = subplot(1, 2, 2);
    surf(X, Y, w_curr); 
    shading interp; colormap(lines(45)); axis off;
    zlim([-1, 1])
    
    
    file_name = [pwd '\images\cool_' num2str(j, '%04d') '.png'];
    saveas(gca, file_name);
    
end

