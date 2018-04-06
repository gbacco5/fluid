% FLUID
% Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath('../../draw','../../tools','../..');


%% DATA
rotor.p = 3; % number of pole pairs
mm = 1e-3; % millimeters
rotor.De = 200*mm; % [m], rotor outer diameter

rotor.Nb = 2; % number of flux-barriers
rotor.tb = [5 10]*mm; % flux-barrier widths
rotor.wc = [6 14 10]*mm; % flux-carrier widths
rotor.Nstep = [2,2]; % number of steps to draw the flux-barrier side
rotor.wrib_t = 1*mm; % [m], tangential iron rib width

% example 3
rotor.barrier_end = 'rect'; % choose 'rect' or comment
rotor.wrib = [1,1]*mm; % [m], radial iron rib widths
rotor.wm = [20,50]*mm;



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
