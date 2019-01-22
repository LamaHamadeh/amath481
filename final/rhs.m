function psift_vec = rhs(t, psif_vec, Lap, A, B, x_vec, y_vec, z_vec, n)

psif = reshape(psif_vec, n, n, n);
psi = ifftn(psif);
psi_vec = reshape(psi, n^3, 1);

spatialPart_vec = ...
    (A(1).*sin(x_vec).^2 + B(1)) .* ...
    (A(2).*sin(y_vec).^2 + B(2)) .* ...
    (A(3).*sin(z_vec).^2 + B(3));

nlPart_vec = (-psi_vec.*conj(psi_vec) + spatialPart_vec).*psi_vec;
nlPart = reshape(nlPart_vec, n, n, n);
nlPartf = fftn(nlPart);

lPartf = (1/2).*Lap.*psif;

psift = 1i*(lPartf + nlPartf);
psift_vec = reshape(psift, n^3, 1);

end
