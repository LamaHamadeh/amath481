function [uvf0_vec] = initialize_uvf(X, Y, m, n)

U0 = tanh(sqrt(X.^2 + Y.^2)).*cos(m*angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2))); 
V0 = tanh(sqrt(X.^2 + Y.^2)).*sin(m*angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));

Uf0 = fft2(U0);
Vf0 = fft2(V0);

uf0_vec = reshape(Uf0, n^2, 1);
vf0_vec = reshape(Vf0, n^2, 1);

uvf0_vec = [uf0_vec; vf0_vec];


end