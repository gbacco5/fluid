% FLUID
% Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath('../../draw','../../tools','../../');


%% DATA
rotor.p = 2; % number of pole pairs
mm = 1e-3; % millimeters
rotor.De = 200*mm; % [m], rotor outer diameter

rotor.Nb = 3; % number of flux-barriers
rotor.tb = [4 8 15]*mm; % flux-barrier widths
rotor.wc = [3 7 12 10]*mm; % flux-carrier widths
rotor.Nstep = [2, 4, 6]; % number of steps to draw the flux-barrier side
rotor.wrib_t = 1*mm; % [m], tangential iron rib width

% example 2
rotor.wrib = [1,1,1]*mm; % [m], radial iron rib widths
rotor.wm = [10,20,40]*mm;



%% barrier points computation
barrier = calc_fluid_barrier(rotor);


%% FEMM drawing
try
  openfemm(1)
  newdocument(0);
  
  draw_fluid_barrier(barrier);
  
catch
  disp('FEMM not available.');
  
end
