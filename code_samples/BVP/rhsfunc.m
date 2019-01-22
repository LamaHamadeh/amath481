%%% RHS function %%%

function f = rhsfunc(t,y,n0,beta)

f1 = y(2);

f2 = (beta-n0)*y(1);

f = [f1;f2];