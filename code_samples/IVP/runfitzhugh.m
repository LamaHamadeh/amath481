% action potential - change 

% a = 0.1; periodic
% a = 0.2; subthreshold
% a = 0.3; single spike
a = 0.3;
b = 0.01;
c = 0.01;

% input current
I=0.1;

vinit = 0.2;
winit =0;


T=400;
deltaT = 0.05;

%set error bounds for integration
options=odeset('RelTol',10^(-10),'AbsTol',[10^(-10)*ones(1,2)]);

% use ode45 to solve the odes
[t,result] = ode15s(@(t,var) fitzhugh(t,var,a,b,c,I),[0:deltaT:T-deltaT],[vinit winit],options);

v = result(:,1);
w = result(:,2);

% refractory var as a function of time
figure(1);
set(gca,'FontSize',18);
hold on; box on;
plot(t,w,'Color',[0 0 1]);
xlabel('t');ylabel('w')

% voltage as a function of time
figure(2);
set(gca,'FontSize',18);
hold on; box on;
plot(t,v,'Color',[0 0 1]);
xlabel('t');ylabel('v')

% phase space
figure(3);
set(gca,'FontSize',18);
plot(v,w);

% draw nulclines
hold on;
vval = [-2:0.01:2];

wl1 = -vval.^3 + (1+a)*vval.^2 -a*vval +I;
wl2 = b/c*vval;

plot(vval,wl1,'r');
plot(vval,wl2,'m');

axis([-0.4 1.2 -0.1 0.4])

xlabel('v');ylabel('w')
