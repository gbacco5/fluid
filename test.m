% TEST
% Free Fluid Flux-Barriers Rotor for Synchronous Reluctance Motor Drawing
%
% Bacco, Giacomo 2018

clear all; close all; clc;
addpath(genpath('./'))
global deb
deb = 0;
tic

if ~exist('test','dir')
  mkdir('test')
end

%% DATA
openfemm(1)
mm = 1e-3; % millimeters
rotor.De = 200*mm; % [m], rotor outer diameter
rotor.wrib_t = 1*mm; % [m], tangential iron rib width
rotor.barrier_end = 'rect'; % choose 'rect' or comment

for p = [1:50] % number of pole pairs
  rotor.p = p;
  for Nb = [1:5] % number of flux-barriers
    
    rotor.Nb = Nb;
    rotor.wrib = zeros(1,Nb);
    % manual input
%     rotor.tb = sqrt(2)/sqrt(rotor.p)*[4 8 15]*mm; % flux-barrier widths
%     rotor.wc = sqrt(2)/sqrt(rotor.p)*[3 7 12 10]*mm; % flux-carrier widths
%     rotor.Nstep = 2*[2, 4, 6]; % number of steps to draw the flux-barrier side    
%     rotor.barrier_angles_el = [14,26,38]*2; % [deg], electrical flux-barrier angles
    
    % 
    if Nb == 1
      rotor.tb = sqrt(2)/sqrt(rotor.p)*[30]*mm; % flux-barrier widths
      rotor.wc = sqrt(2)/sqrt(rotor.p)*[10 22]*mm; % flux-carrier widths
      rotor.Nstep = 2*[6]; % number of steps to draw the flux-barrier side
    elseif Nb == 2
      rotor.tb = sqrt(2)/sqrt(rotor.p)*[10 19]*mm; % flux-barrier widths
      rotor.wc = sqrt(2)/sqrt(rotor.p)*[5 10 18]*mm; % flux-carrier widths
      rotor.Nstep = 2*[3, 5]; % number of steps to draw the flux-barrier side
    elseif Nb == 3
      rotor.tb = sqrt(2)/sqrt(rotor.p)*[4 8 15]*mm; % flux-barrier widths
      rotor.wc = sqrt(2)/sqrt(rotor.p)*[3 7 12 10]*mm; % flux-carrier widths
      rotor.Nstep = 2*[2, 4, 6]; % number of steps to draw the flux-barrier side
    elseif Nb == 4
      rotor.tb = sqrt(2)/sqrt(rotor.p)*[3 6 9 12]*mm; % flux-barrier widths
      rotor.wc = sqrt(2)/sqrt(rotor.p)*[2 4 6 8 10]*mm; % flux-carrier widths
      rotor.Nstep = 2*[3, 4, 6, 8]; % number of steps to draw the flux-barrier side
    elseif Nb == 5
      rotor.tb = sqrt(2)/sqrt(rotor.p)*[2 3 5 6 7]*mm; % flux-barrier widths
      rotor.wc = sqrt(2)/sqrt(rotor.p)*[2 3 4 5 7 8]*mm; % flux-carrier widths
      rotor.Nstep = 2*[3, 4, 6, 8, 10]; % number of steps to draw the flux-barrier side
    end

    
    
    %% barrier points computation
    barrier = calc_fluid_barrier(rotor);
    
        
    %% FEMM drawing
    newdocument(0);
    
    draw_fluid_barrier(barrier)
    mi_saveas(['test/Nb_',num2str(Nb),'_p_',num2str(p),'.fem']);
    mi_close;
    
  end
end

closefemm;
toc
