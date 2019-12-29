% script to solve the Lippmann-Schinwer using the slow formualtion

% discretization and frequency
h = 2.0^(-12);
omega = 1024.0;

% size of box
a = 1;
b = 0.0312;

% discretization
x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

X = repmat(x', 1, m);
Y = repmat(y, n,1);

% the heterogeneity
window = @(y,alpha, beta) 1*(abs(y)<=beta) + ...
                         (abs(y)>beta).*(abs(y)<alpha).* ...
                         exp(2*exp(-(alpha- beta)./(abs(y)-beta))./ ...
                         ((abs(y)-beta)./(alpha- beta)-1) );

nu = @(x,y)  -0.05*(sin(4*pi*x/(0.96))).* ...
                    window(y,0.96*b/2, 0.48*b/2).* ...
                    window(x,0.96*a/2, 0.35);
                     
% the slowness squared
M = 1 + nu(X,Y);

% plotting the slowness squared
% the perturbation oscillated very slowly in the x direction
figure(1); clf();
title("slowness squared")
imagesc(M)

% we define the Lippmann-Schwinger operator
LS = LippmannSchwinger(x,y,omega,nu,a);

% building the right hand-side
u_inc = exp(omega*1i*X);
rhs = -omega^2*nu(X,Y).*u_inc ;

% solving the Lippmann-Schwinger equation
sigma = LS\rhs(:);

% computing the wavefield
u = LS.apply_Green(sigma); 

% plotting the scattered wavefield
figure(2); clf();
imagesc(real(reshape(u, n, m)));

% plotting the total wavefield
figure(3); clf();
imagesc(real(reshape(u+u_inc, n, m)));


% computing the mu solution 
% sigma = mu*exp(omega*1i*(e(1)*X + e(2)*Y ));

LS_Slow = LippmannSchwingerSlow(x,y,omega,nu,a, [1 0]);

rhsSlow = -omega^2*nu(X,Y);

mu = LS_Slow\rhsSlow(:);

% plotting the scattered wavefield
figure(4); clf();
imagesc(abs(reshape(sigma.*conj(u_inc(:))-mu , n, m)) );
colorbar();
title("error between mu slow and sigma/e^{i omega e \cdot x}")

% plotting mu slow
figure(5); clf();
imagesc(real(reshape(mu , n, m)) );
colorbar();
title("mu slow")
