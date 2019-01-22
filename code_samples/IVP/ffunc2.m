function f = ffunc2(t,y,lambda)


f1 = lambda*y(1);
f2 = y(2)-1/2*y(2)^3;

f = [f1;f2];