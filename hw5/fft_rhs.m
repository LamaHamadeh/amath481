function [uvft_vec] = fft_rhs(t, uvf_vec, n, Lap, D1, D2, beta)

uf_vec = uvf_vec(1:n^2);
vf_vec = uvf_vec((n^2 + 1):end);

Uf = reshape(uf_vec, n, n);
Vf = reshape(vf_vec, n, n);

U = real(ifft2(Uf));
V = real(ifft2(Vf));

u_vec = reshape(U, n^2, 1);
v_vec = reshape(V, n^2, 1);

A_vec = u_vec.^2 + v_vec.^2;
lambda_vec = 1 - A_vec;
omega_vec = -beta*A_vec;

utLpart_vec = lambda_vec.*u_vec - omega_vec.*v_vec;
vtLpart_vec = omega_vec.*u_vec + lambda_vec.*v_vec;

utfLpart = fft2(reshape(utLpart_vec, n, n));
vtfLpart = fft2(reshape(vtLpart_vec, n, n));

uft = utfLpart + D1*Lap.*Uf;
vft = vtfLpart + D2*Lap.*Vf;

uvft_vec = [reshape(uft, n^2, 1); reshape(vft, n^2, 1)];

end