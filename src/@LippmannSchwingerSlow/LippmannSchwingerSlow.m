classdef LippmannSchwingerSlow < LippmannSchwinger
    % TO DO: make this work in a rectangle! 
    % this class encodes the fast application of 
    % M = I + omega^2 G * spddiagm(nu)
    % all the properties are public
    properties (SetAccess = public)
        e  % direction of the propagation
        X
        Y
    end
     
    methods
        % constructor
        function LS = LippmannSchwingerSlow(x,y,omega,nu,a,e)
            % x is a vector
            % y is a vector 
            % omega is a real
            % nu is either a function handle or a vector
            
            % suing the constructor of the first class
            LS@LippmannSchwinger(x,y,omega,nu,a)
            LS.e = e;
            
            % this one needs the extra information
            n = length(x); m = length(y);
            LS.X = repmat(x', 1, m);
            LS.Y = repmat(y, n,1);
            
        end
        
        
    end
end