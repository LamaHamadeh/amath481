clear all; close all; % clear previous figures and values
% initialize grid size, time
Time=4; L=20; n=100;
x2=linspace(-L/2,L/2,n+1); 
x=x2(1:n); 
dx=x(2)-x(1)
dt=0.01;
time_steps=Time/dt; 
t=0:dt:Time;

%%%%%
% CFL number 
CFL=dt/(dx^2);
%%%%%

% initial conditions 
u0=exp(-x.^2)'; 
u1=u0+dt*exp(-(x).^2).'; 
usol(:,1)=u0; 
usol(:,2)=u1;

%%%%%
% sparse matrix for second derivative term 
e1=ones(n,1); 
A=spdiags([e1 -2*e1 e1],[-1 0 1],n,n);
A(1,n)=1; A(n,1)=1;
%%%%%

% leap frog (2,2) or euler iteration scheme 
for j=1:time_steps-1 
%%%%%
%---euler---
%u2 = u1 + CFL*A*u1; 
%u1=u2;	
%---euler---
%%%%%

%%%%%
%---leap frog (2,2)---
u2 = u0 + 2*CFL*A*u1;
u0 = u1; u1 = u2;	
%---leap frog (2,2)---
%%%%%

usol(:,j+2)=u2;
end

% plot the data
waterfall(x,t(1:5:end),usol(:,1:5:end)'); map=[0 0 0]; colormap(map);

% set x and y limits and fontsize 
set(gca,'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'FontSize',[20]);
set(gca,'Ylim',[0 Time],'Ytick',[0 Time/2 Time],'FontSize',[20]); 
view(25,40);
% set axis labels and fonts 
xl=xlabel('x'); yl=ylabel('t'); zl=zlabel('u'); 
set(xl,'FontSize',[20]);set(yl,'FontSize',[20]);
set(zl,'FontSize',[20]);

%print -djpeg -r0 fig.jpg % print jpeg at screen resolution