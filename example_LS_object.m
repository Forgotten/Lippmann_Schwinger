% script to solve the Lippmann-Schinwer equation

h = 2.0^(-10);
omega = 256.0;

% size of box
a = 1;
b = 0.25;

% discretization
x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

X = repmat(x', 1, m);
Y = repmat(y, n,1);

% the heterogeneity
nu = @(x,y) -0.3*exp(-320*(x.^2 + y.^2)).*(abs(x)<0.48).*(abs(y)<0.48);

% the slowness squared
M = 1 + nu(X,Y);

% plotting the slowness squared
figure(1); clf();
title("slowness squared")
imagesc(M)

% we define the Lippmann-Schwinger operator
LS = LippmannSchwinger(x,y,omega,nu,a);

% building the right hand-side
u_inc = exp(omega*1i*X);
rhsDual = -omega^2*nu(X,Y).*u_inc ;

% solving the Lippmann-Schwinger equation
sigma = LS\rhsDual(:);

% computing the wavefield
u = LS.apply_Green(sigma); 

% plotting the scattered wavefield
figure(2); clf();
imagesc(real(reshape(u, n, m)));

% plotting the total wavefield
figure(3); clf();
imagesc(real(reshape(u+u_inc, n, m)));