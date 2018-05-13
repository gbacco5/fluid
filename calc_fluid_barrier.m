function barrier = calc_fluid_barrier(r)
% CALC_FLUID_BARRIER computes the flux-barrier points along the streamline
% function.


%% DATA
global deb

Dr = r.De; % [m], rotor outer diameter
ScalingFactor = 1/( 10^(round(log10(Dr))) );
% ScalingFactor = 1;
Dr = Dr*ScalingFactor;

p = r.p; % number of pole pairs
Nb = r.Nb; % number of flux-barriers
tb = r.tb*ScalingFactor; % flux-barrier widths
wc = r.wc*ScalingFactor; % flux-carrier widths
Nstep = r.Nstep; % number of steps to draw the flux-barrier side

wrib_t = r.wrib_t*ScalingFactor; % [m], tangential iron rib width

if isfield(r,'barrier_angles_el')
  barrier_angles_el = r.barrier_angles_el; % [deg], electrical flux-barrier angles
  AutoBarrierEndCalc = 0;
else
  barrier_angles_el = zeros(1,Nb);
  AutoBarrierEndCalc = 1;
end

if isfield(r,'wm')
  wm = r.wm*ScalingFactor;
else
  wm = 0;
end
if isfield(r,'wrib')
  wrib = r.wrib*ScalingFactor + wm; % [m], radial iron rib widths
else
  wrib = zeros(1,Nb) + wm;
end 


Dend = Dr - 2*wrib_t; % [m], flux-barrier end diameter
Dsh = Dend - 2*( sum(tb) + sum(wc) ); % [m], shaft diameter
R0 = Dsh/2; % [m], shaft radius
barrier_angles = barrier_angles_el/p; % [deg], flux-barrier angles
if isfield(r,'barrier_end')
  barrier_end = r.barrier_end;
else
  barrier_end = '';
end


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

vr = @(r,th,R0) p*(r.^(p-1) - R0^(2*p)./r.^(p+1)).*cos(p*th);
vt = @(r,th,R0) -p*(r.^(p-1) + R0^(2*p)./r.^(p+1)).*sin(p*th);
vx = @(vr_v,vth_v,th) vr_v.*cos(th) - vth_v.*sin(th);
vy = @(vr_v,vth_v,th) vr_v.*sin(th) + vth_v.*cos(th);


%% Precomputations
rho0 = rho_map(R0);

%% Central base points
RAprime = Dend(1)/2 - [0, cumsum( tb(1:end-1))] - cumsum(wc(1:end-1)); % top
RBprime = RAprime - tb; % bottom
te_qAxis = pi/(2*p); % q-axis angle in rotor reference frame

% get A' and B' considering rib and magnet widths
mCentral = tan(te_qAxis); % slope
qCentral = repmat( -wrib/2/cos(te_qAxis), 1, 2); % intercept

psiCentralPtA = psi_fluid(rho_map(RAprime), xi_map(te_qAxis), rho0);
psiCentralPtB = psi_fluid(rho_map(RBprime), xi_map(te_qAxis), rho0);
psiCentralPt = [psiCentralPtA, psiCentralPtB];
psiA = psiCentralPtA;
psiB = psiCentralPtB;

CentralPt_Eq = @(th) ...
  r_map( rho_fluid(psiCentralPt, xi_map(th), rho0) ).*...
  ( sin(th) - mCentral*cos(th) ) - qCentral;

if deb == 1
  options.Display = 'iter'; % turn off folve display
else
  options.Display = 'off'; % turn off folve display
end
options.Algorithm = 'levenberg-marquardt'; % non-square systems
options.FunctionTolerance = 1*eps;
options.TolFun = options.FunctionTolerance;
options.StepTolerance = 1e4*eps;
options.TolX = options.StepTolerance;
% I thought the new Matlab syntax for fsolve was options.FunctionTolerance,
% I was wrong.

X0 = repmat(te_qAxis,1,2*Nb);
options = GetFSolveOptions(options);
teAB = fsolve(CentralPt_Eq, X0, options);
teA = teAB(1:Nb);
teB = teAB(Nb+1:end);
RA = r_map( rho_fluid(psiA, xi_map(teA), rho0) );
RB = r_map( rho_fluid(psiB, xi_map(teB), rho0) );

% magnet central base point radius computation
RAsecond = RA.*cos(te_qAxis - teA);
RBsecond = RB.*cos(te_qAxis - teB);  

Rmag = (RAprime + RAsecond + RBprime + RBsecond)/4;


%% Outer base points C,D preparation
RCprime = Dend/2;
teCprime = th_map( xi_fluid(psiA, rho_map(RCprime), rho0) );
xCprime = Dend/2.*cos(teCprime);
yCprime = Dend/2.*sin(teCprime);

RDprime = Dend/2;
teDprime = th_map( xi_fluid(psiB, rho_map(RDprime), rho0) );
xDprime = Dend/2.*cos(teDprime);
yDprime = Dend/2.*sin(teDprime);

if AutoBarrierEndCalc
  teE = (teCprime + teDprime)/2;
  aphE = pi/2/p - teE;
  barrier_angles = 180/pi*aphE;
  barrier_angles_el = p*barrier_angles;
else
  aphE = barrier_angles*pi/180;
  teE = pi/2/p - aphE;
end
xE = Dend/2.*cos(teE);
yE = Dend/2.*sin(teE);  

%% Outer base points C (top)
if strcmp(barrier_end, 'rect')
  RC = RCprime;
  teC = teCprime;
  xC = xCprime;
  yC = yCprime;
  xOC = xC;
  yOC = yC;

else
  options.Algorithm = 'trust-region-dogleg'; % non-square systems
  
  BarrierEndSystem = @(th,xd,yd,xo,yo,R) ...
    [xd - r_map(rho_fluid(psiA', p*th, rho0)).*cos( th )
    yd - r_map(rho_fluid(psiA', p*th, rho0)).*sin( th )
    (xd - xo).^2 + (yd - yo).^2 - R.^2
    (xE' - xo).^2 + (yE' - yo).^2 - R.^2
    (xo - xd).*vx( vr( r_map(rho_fluid(psiA', p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiA', p*th, rho0)) ,th,R0 ), th) + (yo - yd).*vy(  vr( r_map(rho_fluid(psiA', p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiA', p*th, rho0)) ,th,R0 ), th)
    (xo - xE').*yE' - (yo - yE').*xE'
    %   th - xi_fluid((rho_fluid(p*th, psiA', rho0)), psiA', rho0)/p % serve?
    ];
  
  X0 = [ 1.5*teE', 0.9*xE', 0.9*yE', 0.8*xE', 0.8*yE', 0.25*xE'];
  % X0 = [ aph_b, 0, 0, 0, 0, 0];
  X = fsolve( @(x) BarrierEndSystem( x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6) ), X0, options);
  
  xOC = X(:,4)';
  yOC = X(:,5)';
  xC = X(:,2)';
  yC = X(:,3)';
  RC = hypot(xC, yC);
  teC = atan2(yC, xC);
end

%% Outer base points D (bottom)
if strcmp(barrier_end, 'rect')
  RD = RDprime;
  teD = teDprime;
  xD = xDprime;
  yD = yDprime;
  xOD = xD;
  yOD = yD;

else
  options.Algorithm = 'levenberg-marquardt'; % non-square systems
  
  BarrierEndSystem = @(th,xd,yd,xo,yo,R) ...
    [xd - r_map(rho_fluid(psiB', p*th, rho0)).*cos( th )
    yd - r_map(rho_fluid(psiB', p*th, rho0)).*sin( th )
    (xd - xo).^2 + (yd - yo).^2 - R.^2
    (xE' - xo).^2 + (yE' - yo).^2 - R.^2
    (xo - xd).*vx( vr( r_map(rho_fluid(psiB', p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiB', p*th, rho0)) ,th,R0 ), th) + (yo - yd).*vy(  vr( r_map(rho_fluid(psiB', p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiB', p*th, rho0)) ,th,R0 ), th)
    (xo - xE').*yE' - (yo - yE').*xE'
    %   th - xi_fluid((rho_fluid(p*th, psi_d, rho0)), psi_d, rho0)/p % serve?
    ];
  
  X0 = [ 0.8*teE', 0.8*xE', 0.8*yE', xE'*.9, yE'*.9, xE'*.2];
  X = fsolve( @(x) BarrierEndSystem( x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6) ), X0, options);
  
  xOD = X(:,4)';
  yOD = X(:,5)';
  xD = X(:,2)';
  yD = X(:,3)';
  RD = hypot(xD, yD);
  teD = atan2(yD, xD);
end

%% Flux-barrier points
% We already have the potentials of the two flux-barrier sidelines
phiA = phi_fluid( rho_map(RA), xi_map(teA), rho0);
phiB = phi_fluid( rho_map(RB), xi_map(teB), rho0);

phiC = phi_fluid( rho_map(RC), xi_map(teC), rho0);
phiD = phi_fluid( rho_map(RD), xi_map(teD), rho0);


%% Code for single Nstep
% dphiAC = (phiC - phiAprime)./Nstep;
% dphiBD = (phiD - phiBprime)./Nstep;
%
% % we create the matrix of potentials phi needed for points intersections
% PhiAC = phiAprime + cumsum( repmat(dphiAC, Nstep - 1, 1) );
% PhiBD = phiBprime + cumsum( repmat(dphiBD, Nstep - 1, 1) );
%
% PhiAC_vec = reshape(PhiAC, numel(PhiAC), 1);
% PhiBD_vec = reshape(PhiBD, numel(PhiBD), 1);
% PsiAC_vec = reshape( repmat( psiA, Nstep-1, 1), numel(PhiAC), 1 );
% PsiBD_vec = reshape( repmat( psiB, Nstep-1, 1), numel(PhiBD), 1 );
%
% % we find all the barrier points along the streamline
% PsiPhi = @(rho,xi, psi,phi, rho0) ...
%   [psi - psi_fluid(rho, xi, rho0)
%    phi - phi_fluid(rho, xi, rho0)];
%
% X0 = [repmat(rho0*1.1, numel(PhiAC_vec), 1), repmat(pi/4, numel(PhiAC_vec), 1)];
% RhoXi_AC = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiAC_vec, PhiAC_vec, rho0 ), X0, options);
% RhoXi_BD = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiBD_vec, PhiBD_vec, rho0 ), X0, options);
%
% R_AC = reshape( r_map(RhoXi_AC(:,1)), Nstep-1, Nb );
% te_AC = reshape( th_map(RhoXi_AC(:,2)), Nstep-1, Nb );
% R_BD = reshape( r_map(RhoXi_BD(:,1)), Nstep-1, Nb );
% te_BD = reshape( th_map(RhoXi_BD(:,2)), Nstep-1, Nb );


%% Code for different Nsteps
% we find all the barrier points along the streamline
PsiPhi = @(rho,xi, psi,phi, rho0) ...
  [psi - psi_fluid(rho, xi, rho0)
  phi - phi_fluid(rho, xi, rho0)];

% barrier(Nb).R_AC = 0;
% barrier(Nb).R_BD = 0;
% barrier(Nb).te_AC = 0;
% barrier(Nb).te_BD = 0;
barrier(Nb) = struct;

for bkk = 1:Nb
  dphiAC = (phiC(bkk) - phiA(bkk))./Nstep(bkk);
  dphiBD = (phiD(bkk) - phiB(bkk))./Nstep(bkk);
  % we create the matrix of potentials phi needed for points intersections
  PhiAC = phiA(bkk) + cumsum( repmat(dphiAC', Nstep(bkk) - 1, 1) );
  PhiBD = phiB(bkk) + cumsum( repmat(dphiBD', Nstep(bkk) - 1, 1) );
  PsiAC = repmat( psiA(bkk), Nstep(bkk)-1, 1);
  PsiBD = repmat( psiB(bkk), Nstep(bkk)-1, 1);

% 1st try
%   X0 = [repmat(rho0*1.1, numel(PhiAC), 1), repmat(pi/4, numel(PhiAC), 1)];
% 2nd try
%   X0 = [repmat(rho0*1.1, numel(PhiAC), 1), repmat(xi_map(teE(bkk)), numel(PhiAC), 1)];
% 3rd try
  X0 = [linspace(rho0, Dend/2, numel(PhiAC))', linspace(pi/4, xi_map(teE(bkk)), numel(PhiAC))'];
  RhoXi_AC = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiAC, PhiAC, rho0 ), X0, options);
  RhoXi_BD = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiBD, PhiBD, rho0 ), X0, options);
  
  R_AC = r_map(RhoXi_AC(:,1));
  te_AC = th_map(RhoXi_AC(:,2));
  R_BD = r_map(RhoXi_BD(:,1));
  te_BD = th_map(RhoXi_BD(:,2));
  
  if deb
    barrier(bkk).R_AC = R_AC/ScalingFactor;
    barrier(bkk).R_BD = R_BD/ScalingFactor;
    barrier(bkk).te_AC = te_AC;
    barrier(bkk).te_BD = te_BD;
  end
  
  % output of points
%   barrier(bkk).Zeta = [...
  Zeta = [...
    % top side
    xE(bkk) + 1j*yE(bkk)
    xOC(bkk) + 1j*yOC(bkk)
    xC(bkk) + 1j*yC(bkk)
    flipud( R_AC.*exp(1j*te_AC) )
    RA(bkk).*exp(1j*teA(bkk))
    % bottom side
    RB(bkk).*exp(1j*teB(bkk))
    R_BD.*exp(1j*te_BD)
    xD(bkk) + 1j*yD(bkk)
    xOD(bkk) + 1j*yOD(bkk)
    xE(bkk) + 1j*yE(bkk)
    ]/ScalingFactor;
  
  barrier(bkk).X = real(Zeta);
  barrier(bkk).Y = imag(Zeta);
  
  % magnet central base point
  barrier(bkk).Rm = Rmag(bkk)/ScalingFactor;
  
end


%% plot
if deb
  
  % draw the rotor
  figure
  hold on
  tt = linspace(0,pi/p,50);
  plot(R0/ScalingFactor*cos(tt), R0/ScalingFactor*sin(tt), 'k');
  plot(Dr/2/ScalingFactor*cos(tt), Dr/2/ScalingFactor*sin(tt), 'k');
  axis equal
  % plot the flux-barrier central point
  plot(RA/ScalingFactor.*exp(1j*teA), 'rd')
  plot(RB/ScalingFactor.*exp(1j*teB), 'bo')
  
  plot(xE/ScalingFactor, yE/ScalingFactor,'ko')
  
  plot(xOC/ScalingFactor, yOC/ScalingFactor,'go')
  plot(xC/ScalingFactor, yC/ScalingFactor,'ro')
  plot(xOD/ScalingFactor, yOD/ScalingFactor,'co')
  plot(xD/ScalingFactor, yD/ScalingFactor,'bo')
  
  %
  % plot(R_AC.*exp(j*te_AC),'r.-')
  % plot(R_BD.*exp(j*te_BD),'b.-')
  
  for bkk = 1:Nb
    % plot flux-barrier sideline points
    plot(barrier(bkk).R_AC.*exp(1j*barrier(bkk).te_AC),'r.-')
    plot(barrier(bkk).R_BD.*exp(1j*barrier(bkk).te_BD),'b.-')
    
    % plot all the complete flux-barrier
    plot(barrier(bkk).X, barrier(bkk).Y, '.-')
  end
  pause(1e-3)
  
end

end