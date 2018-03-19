% FLUID
% Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath('./tools')


%% DATA
rotor.p = 2; % number of pole pairs
mm = 1e-3; % millimeters
rotor.Nb = 3; % number of flux-barriers
rotor.tb = [4 7 15]*mm; % flux-barrier widths
rotor.wc = [3 6 12 10]*mm; % flux-carrier widths
rotor.Nstep = [2, 4, 6]; % number of steps to draw the flux-barrier side

rotor.De = 200*mm; % [m], rotor outer diameter
rotor.wrib_t = 1*mm; % [m], tangential iron rib width
rotor.barrier_angles_el = [14,26,38]*2; % [deg], electrical flux-barrier angles

rotor.wrib = [1,2,4]*mm; % [m], radial iron rib widths


%% barrier points computation
barrier = calc_fluid_barrier(rotor);


%% simple matlab plot
figure
hold all
axis equal
plot(barrier(1).X, barrier(1).Y, '.-')
plot(barrier(2).X, barrier(2).Y, '.-')
plot(barrier(3).X, barrier(3).Y, '.-')


%% FEMM drawing

openfemm(1)
newdocument(0);

draw_fluid_barrier(barrier)

