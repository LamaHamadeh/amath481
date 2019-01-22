close all;clear variables;

% differentiation cheb matrix +b.c.
N=20;
[D,x] = cheb(N-1);
D2= D^2;
D2(1,:) = zeros(1,N);
D2(N,:) = zeros(1,N);

y=x;
[X,Y]=meshgrid(x,y);

Uinit = exp(-(X.^2+Y.^2)/0.1);
%Uinit = exp(-(X.^2+Y.^2)/0.1).*cos((X-1)*5).*cos((Y-1)*5);
%pcolor(X,Y,Uinit)

%reshape to a vector
uinitvec = reshape(Uinit,N^2,1);
I = eye(length(D2));

% laplacian
Lap = kron(D2,I)+kron(I,D2);

tspan = [0:0.05:0.2];

[t,usolvec] =ode45(@(t,uvec) rhs_heat_cheb(t,uvec,Lap),tspan,uinitvec);


for j=1:length(t)
   
    curU = reshape(usolvec(j,:),N,N);
    pcolor(X,Y,curU);shading interp;
    
    drawnow;
    
    pause(0.1);

    
end