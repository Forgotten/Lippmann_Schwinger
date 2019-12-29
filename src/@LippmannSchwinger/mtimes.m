function LSx = mtimes(LS,x)
% MTIMES  Overload multiplication operator for Ham class
%    HX = H*X returns a wavefun corresponding to H*X.
%
%    See also Ham, Wavefun.

%  Copyright (c) 2019-2020 Leonardo Zepeda-Núñez

% we need to assert that x is a vector and that is has the correct size
assert(size(x,1) == LS.n*LS.m )

BExt = zeros(LS.ne,LS.me);
BExt(1:LS.n,1:LS.m)= reshape(x,LS.n,LS.m);

% Fourier Transform
BFft = fftshift(fft2(BExt));
% Component-wise multiplication
BFft = LS.GFFT.*BFft;
% Inverse Fourier Transform
BExt = ifft2(ifftshift(BFft));

% multiplication by omega^2
B = BExt(1:LS.n, 1:LS.m);

LSx = -x + LS.omega^2*LS.nu.*B(:);

end