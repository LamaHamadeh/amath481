clear all; close all; clc;

beta = 99;
xinit = linspace(-1,1,50);
init = bvpinit(xinit,@bc_init,beta);


sol = bvp4c(@rhs_bvp,@bvp_bc,init);

x = linspace(-1,1,100);
BS = deval(sol,x);
plot(x,BS(1,:));