
function [auto] = get_auto(im)

[Nv,Nu,blank] = size(im);

% calculate 2D autocorrel
im=im(1:min(Nu,Nv),1:min(Nu,Nv));
% 2D-FFT transform on de-meaned image
% power spectrum
mag=abs(fft2(fftshift(im-mean(im(:))))).^2;
%Shift zero-frequency component to centre of spectrum
auto=fftshift(real(ifft2(mag)));
auto = auto./max(auto(:));

[centx,centy] = find(auto==1);
% spectify number of lags to compute
l = length(auto);
nlags=round(l/8);
% centre 2d autocorrelogram
auto = auto(centx-nlags:centx+nlags,centy-nlags:centy+nlags);