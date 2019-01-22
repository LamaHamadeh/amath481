% 1. Create the matrix A

clear all; close all;
n=100;
e1=ones(n,1); % build a vector of ones 
A=spdiags([e1 -2*e1 e1],[-1 0 1],n,n); % diagonals 
A(1,n)=1; A(n,1)=1; % periodic boundaries

% 2. Generate the desired initial condition vector u = u0. 

L=40;
x2=linspace(-L/2,L/2,n+1);

xspan=x2(1:n);	

u0 = exp(-xspan.^2);


% 3. Evolve in time

dx = xspan(2)-xspan(1);
tspan = [0:0.1:10];
k=0.7;

[t,u]=ode45(@(t,u) rhspde(t,u,k,dx,A),tspan,u0);



% 4. Plot the results as a function of time and space.

waterfall(xspan,t,u)

% set colormap and view angle 
map=[0 0 0]; 
colormap(map);
view(25,40)

% set axis labels and fonts 
xl=xlabel('x'); yl=ylabel('t'); zl=zlabel('u','Rotation',0);
set(xl,'FontSize',[20]);set(yl,'FontSize',[20]);set(zl,'FontSize',[20]);

