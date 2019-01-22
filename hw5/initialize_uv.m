function [uv0_vec] = initialize_uv(X, Y, m, n)

U0 = tanh(sqrt(X.^2 + Y.^2)).*cos(m*angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2))); 
V0 = tanh(sqrt(X.^2 + Y.^2)).*sin(m*angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));

u0_vec = reshape(U0, n^2, 1);
v0_vec = reshape(V0, n^2, 1);

uv0_vec = [u0_vec; v0_vec];

end