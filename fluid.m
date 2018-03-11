% FLUID
% Free Fluid Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
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

%% IMPLICIT FUNCTIONS
% definition of fluid past a cylinder functions
psi_fluid = @(rho,xi,rho0) (rho.^2 - rho0^2)./rho.*sin(xi);
phi_fluid = @(rho,xi,rho0) (rho.^2 + rho0^2)./rho.*cos(xi);
xi_fluid  = @(psi,rho,rho0) asin(psi.*rho./(rho.^2 - rho0^2));
rho_fluid = @(psi,xi,rho0) ( psi + sqrt(psi.^2 + 4*sin(xi).^2*rho0^2) )./(2*sin(xi));

%%





%% PLOT
% draw the rotor
figure
hold on
tt = linspace(0,pi/p,50);
plot(R0*cos(tt), R0*sin(tt), 'k');
plot(Dr/2*cos(tt), Dr/2*sin(tt), 'k');
axis equal

RA = Dend/2 - [0, cumsum( tb(1:end-1))] - cumsum(wc(1:end-1));
RB = RA - tb;
plot(RA.*exp(j*pi/2/p), 'rd')
plot(RB.*exp(j*pi/2/p), 'bo')

aph_b = barrier_angles*pi/180;
te_b = pi/2/p - aph_b;
xb = Dend/2*cos(te_b);
yb = Dend/2*sin(te_b);
plot(xb, yb,'ko')











toc