classdef LippmannSchwinger
    % TO DO: make this work in a rectangle! 
    % this class encodes the fast application of 
    % M = I + omega^2 G * spddiagm(nu)
    % all the properties are public
    properties (SetAccess = public)
        
        % array that contains the 
        GFFT
        % array nxm with the compactly supported heterogeneity
        nu
        % number of points in the extended domain
        ne
        me
        % number of points in the physical domain
        n
        m
        % frequency
        omega
        % we can add here a string that encodes the different 
        % quadrature rules
        quadRule
    end
     
    methods
        % constructor
        function LS = LippmannSchwinger(x,y,omega,nu,a)
            % x is a vector
            % y is a vector 
            % omega is a real
            % nu is either a function handle or a vector
            
            % by default we will use Greengard Vico
            
            n = length(x);
            m = length(y);
            
            X = repmat(x', 1, m);
            Y = repmat(y, n,1);
            
            Lp = 4*a ; % we need something at least 4 times bigger
            L  = a*1.5;

            % n needs to be even
            kx = (-(2*n):1:(2*n-1));
            ky = (-(2*m):1:(2*m-1));

            KX = (2*pi/Lp)*repmat(kx', 1, 4*m);
            KY = (2*pi/Lp)*repmat(ky, 4*n,1);

            S = sqrt(KX.^2 + KY.^2);

            G2D = @(L,k,s)(1 + ...
                           (1i*pi/2*L*besselh(0,1,L*k)).*(s.*besselj(1,L*s)) - ...
                           (1i*pi/2*L*k*besselh(1,1,L*k)).*besselj(0,L*s)...
                           )./(s.^2 - k^2);

            LS.GFFT = G2D(L, omega, S);

% degrees of freedom of the extended fft

            LS.n = n;
            LS.m = m;
            LS.ne = 4*n;
            LS.me = 4*m;
            LS.omega = omega;
            
            nu_vect = nu(X,Y);
            LS.nu = nu_vect(:);
            % in some cases we 
            if ~isempty(isnan(LS.nu))
                fprintf("some values of the window may be Nan")
                LS.nu(isnan(LS.nu)) = 0.0;
            end
        end
        
        % we need to use this one for the multiplication
        function Gu = apply_Green(LS, u)
            
            BExt = zeros(LS.ne,LS.me);
            BExt(1:LS.n,1:LS.m)= reshape(u,LS.n,LS.m);

            % Fourier Transform
            BFft = fftshift(fft2(BExt));
            % Component-wise multiplication
            BFft = LS.GFFT.*BFft;
            % Inverse Fourier Transform
            BExt = ifft2(ifftshift(BFft));

            % multiplication by omega^2
            Gu = BExt(1:LS.n, 1:LS.m);
        end
        
    end
end