I = imread('1noise.jpg');

[H W D] =size(I);

figure;imshow(I);



for j=1:3
F = fft2(I(:,:,j));
F= fftshift(F);

%imshow(abs(F),[0,10^5]);
F1 = zeros(H,W);
F1(15:end-14,25:end-24) = F(15:end-14,25:end-24);

F1 = ifftshift(F1);
F1 = ifft2(F1);


I(:,:,j) = imadjust(uint8(abs(F1)));

end

%I = histeq(I);

figure;imshow(I);

