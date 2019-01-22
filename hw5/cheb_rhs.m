function [uvt_vec] = cheb_rhs(t, uv_vec, n, Lap, D1, D2, beta)

u_vec = uv_vec(1:n^2);
v_vec = uv_vec((n^2 + 1):end);

A_vec = u_vec.^2 + v_vec.^2;

lambda_vec = 1 - A_vec;
omega_vec = -beta*A_vec;

ut_vec = lambda_vec.*u_vec - omega_vec.*v_vec + D1*Lap*u_vec;
vt_vec = omega_vec.*u_vec + lambda_vec.*v_vec + D2*Lap*v_vec;

uvt_vec = [ut_vec; vt_vec];

end