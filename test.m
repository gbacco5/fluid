% TEST
% Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath(genpath('./'))
global deb
deb = 1;

%% DATA
openfemm(1)
mm = 1e-3; % millimeters
rotor.De = 200*mm; % [m], rotor outer diameter
rotor.wrib_t = 1*mm; % [m], tangential iron rib width
% rotor.barrier_end = 'rect'; % choose 'rect' or comment
rotor.wrib = 0*[1,2,4]*mm; % [m], radial iron rib widths
 rotor.barrier_end = 'rect'; % choose 'rect' or comment

for p = 4 % number of pole pairs
  rotor.p = p;
  for Nb = 3 % number of flux-barriers
    
    rotor.Nb = Nb;
     rotor.tb = sqrt(2)/sqrt(rotor.p)*[4 8 15]*mm; % flux-barrier widths
     rotor.wc = sqrt(2)/sqrt(rotor.p)*[3 7 12 10]*mm; % flux-carrier widths
%     rotor.tb = 2/(rotor.p)*[4 8 15]*mm; % flux-barrier widths
%     rotor.wc = 2/(rotor.p)*[3 7 12 10]*mm; % flux-carrier widths
%    rotor.tb = [2 3 4]*mm; % flux-barrier widths
%    rotor.wc = [2 3 4 5]*mm; % flux-carrier widths
    rotor.Nstep = 2*[2, 4, 6]; % number of steps to draw the flux-barrier side
    
    rotor.barrier_angles_el = [14,26,38]*2; % [deg], electrical flux-barrier angles
    
    
    
    %% barrier points computation
    barrier = calc_fluid_barrier(rotor);
    
        
    %% FEMM drawing
    
    newdocument(0);
    
    draw_fluid_barrier(barrier)
    mi_saveas(['test/p_',num2str(p),'.fem']);
    mi_close;
    
  end
end
