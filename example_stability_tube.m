% starting up the path
LS_startup()
% discretization and frequency
h = 2.0^(-12);
nppw = 4 ;
omega = (1/h)/nppw;

% size of box
a = 1;
b = 1;

% discretization
x = -a/2:h:a/2-h;
y = -b/2:h:b/2-h;

n = length(x); m = length(y);

X = repmat(x', 1, m);
Y = repmat(y, n,1);

% the window
window = @(y,alpha, beta)1*(abs(y)<=beta) + (abs(y)>beta).*(abs(y)<alpha).* ...
                         exp(2*exp(-(alpha-beta)./(abs(y)-beta))./ ...
                         ((abs(y)-beta)./(alpha-beta)-1) );
% the heterogeneity 
mm = @(x,y)  -0.05*(sin(4*pi*x/(0.96))).* ...
                    window(y,0.960001*b/2, 0.480001*b/2).* ...
                    window(x,0.960001*a/2, 0.350001);
                     
% the slowness squared
nsqr = 1 + mm(X,Y);

% plotting the slowness squared
% the perturbation oscillated very slowly in the x direction
figure(1); clf();
title("slowness squared")
imagesc(nsqr)

% we define the Lippmann-Schwinger operator
LS = LippmannSchwinger(x,y,omega,mm,a);

% building the right hand-side
u_inc = exp(omega*1i*X);
m_vect = mm(X,Y);
m_vect(isnan(m_vect)) = 0.0;
rhs = -omega^2*m_vect.*u_inc ;

% solving the Lippmann-Schwinger equation
sigma = LS\rhs(:);

% computing the wavefield
u = LS.apply_Green(sigma); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% solving the same problem but with a much thinner tube
% size of box
a1 = 1;
b1 = 0.0078125;

% discretization
x1 = -a1/2:h:a1/2-h;
y1 = -b1/2:h:b1/2-h;

n1 = length(x1); m1 = length(y1);

X1 = repmat(x1', 1, m1);
Y1 = repmat(y1, n1,1);

% the window
window = @(y,alpha, beta)1*(abs(y)<=beta) + (abs(y)>beta).*(abs(y)<alpha).* ...
                         exp(2*exp(-(alpha-beta)./(abs(y)-beta))./ ...
                         ((abs(y)-beta)./(alpha-beta)-1) );
% the heterogeneity 
mm = @(x,y)  -0.05*(sin(4*pi*x/(0.96))).* ...
                    window(y,0.960001*b1/2, 0.480001*b1/2).* ...
                    window(x,0.960001*a1/2, 0.350001);
                     
% the slowness squared
nsqr = 1 + mm(X1,Y1);

% plotting the slowness squared
% the perturbation oscillated very slowly in the x direction
figure(2); clf();
title("slowness squared")
imagesc(nsqr)

% we define the Lippmann-Schwinger operator
LS_1 = LippmannSchwinger(x1,y1,omega,mm,a1);

% building the right hand-side
u_inc = exp(omega*1i*X1);
m_vect = mm(X1,Y1);
m_vect(isnan(m_vect)) = 0.0;
rhs = -omega^2*m_vect.*u_inc ;

% solving the Lippmann-Schwinger equation
sigma1 = LS_1\rhs(:);

% computing the wavefield
u1 = LS_1.apply_Green(sigma1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting 

% plotting the scattered wavefield
figure(3); clf();
imagesc(real(reshape(u1, n1, m1)));

% plotting the scattered wavefield
figure(4); clf();
imagesc(real(reshape(u, n, m)));

u_vec = reshape(u, n,m); 
u_vec1 = reshape(u1, n1,m1); 

% plotting the total wavefield
figure(5); clf();
semilogy(abs(u_vec(:,end/2)-u_vec1(:,end/2))/norm(u_vec(:,end/2)),...
'LineWidth', 3);
set(gca,'FontSize',18)
title("error u following the principal ray")

sigma_vec = reshape(sigma, n,m); 
sigma_vec1 = reshape(sigma1, n1,m1); 

% plotting the total wavefield
figure(6); clf();
semilogy(abs(sigma_vec(:,end/2)-sigma_vec1(:,end/2))/norm(sigma_vec(:,end/2)),...
'LineWidth', 3);
set(gca,'FontSize',18) 
title("error sigma following the principal ray")
