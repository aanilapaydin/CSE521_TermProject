function [errout,xo,to,Uo] = heatBTCS(nt,nx,alpha,L,tmax,errPlots)
% heatBTCS Solve 1D heat equation with the BTCS scheme
%
% Synopsis: heatBTCS
% heatBTCS(nt)
% heatBTCS(nt,nx)
% heatBTCS(nt,nx,alpha)
% heatBTCS(nt,nx,alpha,L)
% heatBTCS(nt,nx,alpha,L,tmax)
% heatBTCS(nt,nx,alpha,L,tmax,errPlots)
% err = heatBTCS(...)
% [err,x,t,U] = heatBTCS(...)
%
% Input: nt = number of steps. Default: nt = 10;
% nx = number of mesh points in x direction. Default: nx=20
% alpha = diffusion coefficient. Default: alpha = 0.1
% L = length of the domain. Default: L = 1;
% tmax = maximum time for the simulation. Default: tmax = 0.5
% errPlots = flag (1 or 0) to control whether plots of the solution
% and the error at the last time step are created.
% Default: errPlots = 1 (make the plots)
%
% Output: err = L2 norm of error evaluated at the spatial nodes on last time step
% x = location of finite difference nodes
% t = values of time at which solution is obtained (time nodes)
% U = matrix of solutions: U(:,j) is U(x) at t = t(j)
if nargin<1, nt = 10; end
if nargin<2, nx = 20; end
if nargin<3, alpha = 0.1; end % diffusivity
if nargin<4, L = 1; end
if nargin<5, tmax = 0.5; end
if nargin<6, errPlots=1; end % flag: no error plots if errPlots=0
% --- Compute mesh spacing and time step
dx = L/(nx-1);
dt = tmax/(nt-1);
% --- Create arrays to save data for export
x = linspace(0,L,nx)’;
t = linspace(0,tmax,nt);
U = zeros(nx,nt);
% --- Set IC and BC
U(:,1) = sin(pi*x/L); % implies u0 = 0; uL = 0;
u0 = 0; uL = 0;
% --- Coefficients of the tridiagonal system
a = (-alpha/dx^2)*ones(nx,1); % subdiagonal a: coefficients of phi(i-1)
c = a; % superdiagonal c: coefficients of phi(i+1)
b = (1/dt)*ones(nx,1) - 2*a; % diagonal b: coefficients of phi(i)
b(1) = 1; c(1) = 0; % Fix coefficients of boundary nodes
b(end) = 1; a(end) = 0;
[e,f] = tridiagLU(a,b,c); % Get LU factorization of coefficient matrix

% O(n) * O (tridiagLUSolve of complexity) = > O(n ^ 2)

% --- Loop over time steps
for m=2:nt
	d = U(:,m-1)/dt; % update right hand side
	d(1) = u0; d(end) = uL; % overwrite BC values
	U(:,m) = tridiagLUSolve(d,a,e,f,U(:,m-1)); % solve the system
end


% --- Compare with exact solution at end of simulation
% ... Remainder of code is the same as heatFTCS