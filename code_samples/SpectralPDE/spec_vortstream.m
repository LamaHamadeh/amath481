clear all;close all;clc

% initial setup
tspan = [0:0.5:25];
nu = 0.001;
L=20;n=64;
x2 = linspace(-L/2,L/2,n+1);x=x2(1:n);
y=x;

kx=(2*pi/L)*[0:(n/2-1) (-n/2):-1];kx(1)=10^-6;
ky=kx;

[X,Y]=meshgrid(x,y);
[KX,KY]=meshgrid(kx,ky);

K = KX.^2+KY.^2;
Kvec = reshape(K,n^2,1);

% initial conditions
%winit=exp(-X.^2-(Y.^2/20));
%winit = exp(-0.25*(X).^2-2*Y.^2);
winit = -exp(-0.25*(X-2).^2-2*Y.^2)+exp(-0.25*(X+2).^2-2*Y.^2);

%winit = -exp(-0.25*(X-2).^2-2*(Y+6).^2)...
%        +exp(-0.25*(X+2).^2-2*(Y+6).^2)...
%        +exp(-0.25*(X-2).^2-2*(Y-6).^2)...
%        -exp(-0.25*(X+2).^2-2*(Y-6).^2);
    
% initial conditions in fourier space, vector form
wfinit = fft2(winit);
wfvecinit = reshape(wfinit,n^2,1);

% integrate in fourier space
[t,wfvecsol]=ode45(@(t,wfvec) vortrhs(t,wfvec,nu,K,Kvec,n,KX,KY),tspan,wfvecinit);

% return to solution space and show as movie
for j=1:length(t)
    
    curw=real(ifft2(reshape(wfvecsol(j,:),n,n)));
    pcolor(X,Y,curw);shading interp;
    drawnow;
    
    pause(0.1);
end

%pcolor(X,Y,w);shading interp;