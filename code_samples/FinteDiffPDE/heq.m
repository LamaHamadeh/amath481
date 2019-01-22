clear all; close all; 

% initialize grid size, time
Time=1; L=20; n=100;
x2=linspace(-L/2,L/2,n+1); 
x=x2(1:n); 
dx=x(2)-x(1); 
dt=0.005;
time_steps=Time/dt; 
t=0:dt:Time;

% CFL
CFL = dt/dx^2

usol = zeros(length(x),length(t));

% ic
u0 = exp(-x.^2).';
u1 = exp(-(x+dt).^2).';

usol(:,1) = u0;
usol(:,2) = u1;

% sparse matrix for second derivative term 
e1=ones(n,1); 
A=spdiags([e1 -2*e1 e1],[-1 0 1],n,n);
A(1,n)=1; A(n,1)=1;


for j=2:time_steps-1
    
    %u2 = u1+CFL*A*u1;
    
    %u1=u2;
    
    %usol(:,j+1) = u2;
    
    u2= u0+2*CFL*A*u1;
    u0=u1;u1=u2;
    
    usol(:,j+1)=u2;
    
    
end

waterfall(x,t(1:50:end),usol(:,1:50:end).'); map=[0 0 0]; colormap(map);
    










