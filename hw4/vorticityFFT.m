function rhs = vorticityFFT(t, w_vec, rhsFun, K, n) 

wf = fft2(reshape(w_vec, n, n));

psif = -wf./K;
psi = real(ifft2(psif));
psi_vec = reshape(psi, n^2, 1);

rhs = rhsFun(psi_vec, w_vec);

end