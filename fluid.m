% FLUID
% Free Fluid Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath('./tools')
tic

%% DATA
p = 2; % number of pole pairs
mm = 1e-3; % millimeters
Nb = 3; % number of flux-barriers
tb = [4 7 15]*mm; % flux-barrier widths
wc = [3 6 12 10]*mm; % flux-carrier widths
Nstep = 5; % number of steps to draw the flux-barrier side

Dr = 200*mm; % [m], rotor outer diameter
wrib_t = 1*mm; % [m], tangential iron rib width
Dend = Dr - 2*wrib_t; % [m], flux-barrier end diameter
Dsh = Dend - 2*( sum(tb) + sum(wc) ); % [m], shaft diameter
R0 = Dsh/2; % [m], shaft radius
barrier_angles_el = [14,26,38]*2; % [deg], electrical flux-barrier angles
barrier_angles = barrier_angles_el/p; % [deg], flux-barrier angles

% wm = []
wrib = [1,2,3]*mm; % [m], radial iron rib widths


%% IMPLICIT FUNCTIONS
% definition of fluid past a cylinder functions
psi_fluid = @(rho,xi,rho0) (rho.^2 - rho0^2)./rho.*sin(xi);
phi_fluid = @(rho,xi,rho0) (rho.^2 + rho0^2)./rho.*cos(xi);
xi_fluid  = @(psi,rho,rho0) asin(psi.*rho./(rho.^2 - rho0^2));
rho_fluid = @(psi,xi,rho0) ( psi + sqrt(psi.^2 + 4*sin(xi).^2*rho0^2) )./(2*sin(xi));

r_map = @(rho) rho.^(1/p);
th_map = @(xi) xi./p;
rho_map = @(r) r.^p;
xi_map = @(th) th.*p;


%% Precomputations
rho0 = rho_map(R0);

%% Central base points
RA = Dend/2 - [0, cumsum( tb(1:end-1))] - cumsum(wc(1:end-1)); % top
RB = RA - tb; % bottom
te_qAxis = pi/(2*p); % q-axis angle in rotor reference frame

% get A' and B' considering rib and magnet widths
mCentral = tan(te_qAxis); % slope
qCentral = -wrib/2/cos(te_qAxis); % intercept

psiCentralPtA = psi_fluid(rho_map(RA), xi_map(te_qAxis), rho0);
psiCentralPtB = psi_fluid(rho_map(RB), xi_map(te_qAxis), rho0);

CentralPtA_Eq = @(th) ...
  r_map( rho_fluid(psiCentralPtA, th, rho0) ).*...
  ( sin(th) - mCentral*cos(th) ) - qCentral;

%if isOctave
%  options = [];
%else
  options.Display = 'off'; % turn off folve display
  options.Algorithm = 'levenberg-marquardt'; % non-square systems
  options.FunctionTolerance = 10*eps;
  options.StepTolerance = 1e4*eps;
%end

X0 = [te_qAxis];
%pars.R0 = R0;
%pars.rho0 = rho0;
%pars.p = p;
%pars.mCentral = mCentral;
%pars.qCentral = qCentral;
%pars.psiCentralPtA = psiCentralPtA;
%pars.psi_fluid = psi_fluid;
%pars.rho_fluid = rho_fluid;
%pars.xi_fluid = xi_fluid;
%pars.r_map = r_map;
%pars.th_map = th_map;
%pars.rho_map = rho_map;
%pars.xi_map = xi_map;

% 
X = SolveSystemMatlabOctave(CentralPtA_Eq, X0, options);
% X1 = fsolve(CentralPtA_Eq, X0, options);





% plot
% % draw the rotor
% figure
% hold on
% tt = linspace(0,pi/p,50);
% plot(R0*cos(tt), R0*sin(tt), 'k');
% plot(Dr/2*cos(tt), Dr/2*sin(tt), 'k');
% axis equal
% 
% plot(RA.*exp(j*pi/2/p), 'rd')
% plot(RB.*exp(j*pi/2/p), 'bo')
% 
% aph_b = barrier_angles*pi/180;
% te_b = pi/2/p - aph_b;
% xb = Dend/2*cos(te_b);
% yb = Dend/2*sin(te_b);
% plot(xb, yb,'ko')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
toc