%% AMATH 481 - PS 0 - Spencer Pease
%

clear all; close all;


%% Exercise 1 - Building a Matrix
%

A = [34 45; 17 6];

save A1.dat A -ascii;

%% Exercise 2 - Matrix Operations
%

A = [1 2; -1 1];
B = [2 0; 0 2];
C = [2 0 -3; 0 0 -1];
D = [1 2; 2 3; -1 0];

x = [1; 0];
y = [0; 1];
z = [1; 2; -1];

A2 = A + B;
A3 = 3*x - 4*y;
A4 = A*x;
A5 = B*(x - y);
A6 = D*x;
A7 = D*y + z;
A8 = A*B;
A9 = B*C;
A10 = C*D;

save A2.dat A2 -ascii
save A3.dat A3 -ascii
save A4.dat A4 -ascii
save A5.dat A5 -ascii
save A6.dat A6 -ascii
save A7.dat A7 -ascii
save A8.dat A8 -ascii
save A9.dat A9 -ascii
save A10.dat A10 -ascii


%% Exercise 3 - Root Finding
%

f = @(x) -1*x - cos(x);
df = @(x) -1 + sin(x);

tol = 1e-6;
maxIter = 100;

x1 = -3;
a = -3;
b = 1;

[x_nr, iter_nr] = newtonRaphson(f, df, x1, tol, maxIter);
[x_bi, iter_bi] = bisection(f, a, b, tol, maxIter);

A11 = x_nr;
A12 = x_bi;
A13 = [iter_nr iter_bi];

save A11.dat A11 -ascii
save A12.dat A12 -ascii
save A13.dat A13 -ascii









