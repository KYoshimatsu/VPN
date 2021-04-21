% ``Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'' [J. Comput. Phys. 390 (2019) 452-469]
% and its Corrigendum.
% by T. Sakurai, K. Yoshimatsu, N. Okamoto and K. Schneider
% This is the Matlab/Octave code for Fig. C19 in 
%  "Corrigendum to `Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'
%-----
% This code is based on the Octave/Matlab code for the example shown in section 2.3 of
% "Analysis and discretization of the volume penalized Laplace operator 
% with Neumann boundary conditions" by D. Kolomenskiy, R. Nguyen van yen
% and K. Schneider
clear all;
close all;

% appendix C
alp=1.0;

% Number of grid points
N = 2048;

% Penalization parameter
eta = 1e-8;

% Domain size
Lx = 2*pi;

% Grid step
h = Lx/N;

% Grid
x = (0:N-1).'*h;

% Fluid domain array index
Omf = find((x>0)&(x<pi));
%Omf1 =find((ge(x,0))&(le(x,pi)));
% Right-hand side of the equation

% Mask function
%chi = zeros(N,1);
chiN = zeros(N,1);
chiD = zeros(N,1);
chiN((x>pi)&(x<1.5*pi)) = 1;
chiD(x>1.5*pi) = 1;
% eps(1) is the machine epsilon near 1. 
% Note that there has to be a grid point on the boundary.
%
chiN(abs(x-pi)<2*eps(1)) = 1/2;
chiN(abs(x-1.5*pi)<2*eps(1)) = 1/2;
chiD(abs(x)<2*eps(1)) = 1/2;
chiD(abs(x-1.5*pi)<2*eps(1)) = 1/2;
%
chi=chiN+chiD;

%f = m^2*cos(m*x);
f = (1-chi).*cos(x);

% Exact solution 
w = cos(x)+alp *x-1;
% Its extension in the solid domain
w(x>pi) = 0;

% Discrete operator matrix
Db = toeplitz([1 -1 zeros(1,N-2)].',[1 zeros(1,N-2) -1])/h;
Df = toeplitz([-1 zeros(1,N-2) 1].',[-1 1 zeros(1,N-2)])/h;
Th = diag(1-chiN+eta*chiN);
%Th = diag(1-chi+eta*chi);
A = -1/2*(Df*Th*Db+Db*Th*Df);

% inhomogeneous Neumann b.c.
bet = alp*sin(x)+ones(N,1)*alp;
% beta * dx chi
dTh = 1/2*bet.*(Db+Df)*chiN;
% Dirichlet pena
Dpena= -chiD/eta ;
%
A=A - diag(Dpena);
f = f + dTh;

% Numerical solution
u = A\f;

% Plot the exact and the numerical solution
figure(1);
plot(x,w,'b-',x,u,'r-');

% Compute the error norm
%err_linfty = max(abs(u(Omf)-w(Omf)))
%err_l2 = sqrt(h*sum((u(Omf)-w(Omf)).^2))

err_linfty = max(abs(u(Omf)-w(Omf)));
err_l1 = h*sum(abs(u(Omf)-w(Omf)));
err_l2 = sqrt(h*sum((u(Omf)-w(Omf)).^2));

error(1:5)=[N h err_linfty err_l1 err_l2];
save('err_DN.data','error');
