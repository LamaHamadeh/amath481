function dydt = fitzhugh(t,y,a,b,c,I)

% voltage variable
v = y(1);

% refractory variable
w = y(2);


% define the 2d system
dv = -v^3 + (1+a)*v^2 -a*v -w + I;
dw = b*v-c*w;

dydt = [dv;dw];