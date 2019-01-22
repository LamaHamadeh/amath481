function rhs = vortrhs(t,wfvec,nu,K,Kvec,n,KX,KY)


wf = reshape(wfvec,n,n);

psif = -wf./K;

psi_x = real(ifft2(i*KX.*psif));
psi_y = real(ifft2(i*KY.*psif));
w_x = real(ifft2(i*KX.*wf));
w_y = real(ifft2(i*KY.*wf));



rhs = -nu*Kvec.*wfvec + reshape(fft2(-psi_x.*w_y+psi_y.*w_x),n^2,1);