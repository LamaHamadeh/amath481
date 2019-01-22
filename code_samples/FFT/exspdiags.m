clear all; close all;
n=25;
e1=ones(n,1); % build a vector of ones 

A=spdiags([e1 -2*e1 e1],[-1 0 1],n,n); % diagonals 

A(1,n)=1; A(n,1)=1; % periodic boundaries

spy(A) % spy the matrix