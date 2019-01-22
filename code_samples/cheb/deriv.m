

x = -1:0.01:1;

u  =  exp(x).*sin(5*x);
ux =  exp(x).*sin(5*x) + 5*exp(x).*cos(5*x);
uxx = - 24*sin(5*x).*exp(x) + 10*cos(5*x).*exp(x);

N=20;

[D,x2] = cheb(N);

D2 = D^2;

u2 = exp(x2).*sin(5*x2);
u2x = D*u2;

u2xx = D2*u2;


subplot(3,1,1);
plot(x,u,'k-',x2,u2,'mo');
subplot(3,1,2);
plot(x,ux,'k-',x2,u2x,'mo');
subplot(3,1,3);
plot(x,uxx,'k-',x2,u2xx,'mo');