function LSmu = mtimes(LS,mu)
% MTIMES  Overload multiplication operator for LippmannSchwingerSlow class
%    LSmu = LS*mu returns a vector 
%
%    See also LippmannSchwinger

%  Copyright (c) 2019-2020 Leonardo Zepeda-Núñez

% we need to assert that x is a vector and that is has the correct size
assert(size(mu,1) == LS.n*LS.m )
% it only works for vectors so far
assert(size(mu,2) == 1 )

u_inc = reshape(exp(1i*LS.omega*(LS.e(1)*LS.X+ LS.e(2)*LS.Y)), LS.n*LS.m, 1) ;
b = mu.*u_inc;

Gu = LS.apply_Green(b);

LSmu = -mu + LS.omega^2*(LS.nu.*(conj(u_inc).*Gu(:)));

end