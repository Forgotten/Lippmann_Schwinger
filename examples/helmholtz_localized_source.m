% script to solve the Lippmann-Schinwer equation

h = 2.0^(-10);
omega = 256.0;

% size of box
a = 1;
b = 1;

% discretization
x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

X = repmat(x', 1, m);
Y = repmat(y, n,1);

theta = [1, 0];

% building the perturbation 
nu = @(x,y) exp(-omega*(x.^2 + y.^2)).*(abs(x)<0.48).*(abs(y)<0.48).*...
    exp(1i*omega*(theta(1)*x + theta(2)*y)) ;

% the right-hand side
rhs = nu(X,Y);

% plotting the slowness squared
figure(1); clf();
title("source")
subplot(1,2,1);
imagesc(real(rhs))
subplot(1,2,2);
imagesc(imag(rhs))


% we define the Lippmann-Schwinger operator
LS = LippmannSchwinger(x,y,omega,nu,a);

% we applyt the Green's fucntion of the background 
u = LS.apply_Green(rhs); 

% plotting the scattered wavefield
figure(2); clf();
title("solution")
subplot(1,3,1);
imagesc(real(reshape(u, n, m)));
subplot(1,3,2);
imagesc(imag(reshape(u, n, m)));
subplot(1,3,3);
imagesc(abs(reshape(u, n, m)));

