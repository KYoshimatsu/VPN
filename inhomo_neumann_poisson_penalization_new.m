% ``Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'' [J. Comput. Phys. 390 (2019) 452-469]
% and its Corrigendum.
% by T. Sakurai, K. Yoshimatsu, N. Okamoto and K. Schneider
% This is the Matlab/Octave code for Fig. 8 in 
%  "Corrigendum to `Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry'
%-----
% This code is based on the Matlab code for the example shown in section 2.3 of
% "Analysis and discretization of the volume penalized Laplace operator 
% with Neumann boundary conditions" by D. Kolomenskiy, R. Nguyen van yen
% and K. Schneider

clear all;
close all;

% This script only works if m is odd
% For m is even, the a additive constant is non-zero
m = 1;
alp=0.0;

% Number of grid points
N = 512;

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

% Mask function
chi = zeros(N,1);
chi(x>pi) = 1;
chi(abs(x-pi)<2*eps(1)) = 1/2;
chi(abs(x)<2*eps(1)) = 1/2;
% eps(1) is the machine epsilon near 1. 
% Note that there has to be a grid point on the boundary.

% Right-hand side of the equation
%f = m^2*cos(m*x);
af=m^2*sin(m*x);
f = (1-chi).*af;

% The r.h.s. of the equation ibecomes discontinuous
% when we extend it in the solid domain. Its value at the point of
% discontinuity must be the half-sum of the limits from the right
% and from the left, for consistency with the central discretization
% of the differential operator.

% Exact solution 
%w = cos(m*x);
 w = sin(m*x)+alp*x-pi*alp/2.0-2.0/(pi*m);
% Its extension in the solid domain
w(x>pi) = 0;

% Discrete operator matrix
Db = toeplitz([1 -1 zeros(1,N-2)].',[1 zeros(1,N-2) -1])/h;
Df = toeplitz([-1 zeros(1,N-2) 1].',[-1 1 zeros(1,N-2)])/h;
Th = diag(1-chi+eta*chi);
A = -1/2*(Df*Th*Db+Db*Th*Df);

% inhomogeneous Neumann b.c.
bet = m*cos(m*x)+ones(N,1)*alp;
% beta * dx chi
dTh = 1/2*bet.*(Db+Df)*chi;

f = f + dTh;

% Ensure zero mean

%f(1) = 0;
%A(1,:) = 1;

% add
%A(1,:)=0;
%A(1,Omf) = 1;

% 1 -> 3
f(3) = 0;
A(3, :) = 1-chi;
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
save('err_1Dinhomo.data','error');
