function barrier = calc_fluid_barrier(r)
% CALC_FLUID_BARRIER computes the flux-barrier points along the streamline
% function.


%% DATA
global deb

p = r.p; % number of pole pairs
Nb = r.Nb; % number of flux-barriers
tb = r.tb; % flux-barrier widths
wc = r.wc; % flux-carrier widths
Nstep = r.Nstep; % number of steps to draw the flux-barrier side

Dr = r.De; % [m], rotor outer diameter
wrib_t = r.wrib_t; % [m], tangential iron rib width
barrier_angles_el = r.barrier_angles_el; % [deg], electrical flux-barrier angles
if isfield(r,'wm')
  wm = r.wm;
else
  wm = 0;
end
wrib = r.wrib + wm; % [m], radial iron rib widths


Dend = Dr - 2*wrib_t; % [m], flux-barrier end diameter
Dsh = Dend - 2*( sum(tb) + sum(wc) ); % [m], shaft diameter
R0 = Dsh/2; % [m], shaft radius
barrier_angles = barrier_angles_el/p; % [deg], flux-barrier angles


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
RA = Dend/2 - [0, cumsum( tb(1:end-1))] - cumsum(wc(1:end-1)); % top
RB = RA - tb; % bottom
te_qAxis = pi/(2*p); % q-axis angle in rotor reference frame

% get A' and B' considering rib and magnet widths
mCentral = tan(te_qAxis); % slope
qCentral = repmat( -wrib/2/cos(te_qAxis), 1, 2); % intercept

psiCentralPtA = psi_fluid(rho_map(RA), xi_map(te_qAxis), rho0);
psiCentralPtB = psi_fluid(rho_map(RB), xi_map(te_qAxis), rho0);
psiCentralPt = [psiCentralPtA, psiCentralPtB];
psiA = psiCentralPtA;
psiB = psiCentralPtB;

CentralPt_Eq = @(th) ...
  r_map( rho_fluid(psiCentralPt, th, rho0) ).*...
  ( sin(th) - mCentral*cos(th) ) - qCentral;

options.Display = 'off'; % turn off folve display
options.Algorithm = 'levenberg-marquardt'; % non-square systems
options.FunctionTolerance = 10*eps;
options.StepTolerance = 1e4*eps;

X0 = repmat(te_qAxis,1,2*Nb);
options = GetFSolveOptions(options);
teABprime = fsolve(CentralPt_Eq, X0, options);
teAprime = teABprime(1:Nb);
teBprime = teABprime(Nb+1:end);
RAprime = r_map( rho_fluid(psiA, xi_map(teAprime), rho0) );
RBprime = r_map( rho_fluid(psiB, xi_map(teBprime), rho0) );


%% Outer base points C (top)
aphE = barrier_angles*pi/180;
teE = pi/2/p - aphE;
xE = Dend/2*cos(teE);
yE = Dend/2*sin(teE);

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

X0 = [ teE', 0*xE', 0*yE', 0*xE', 0*yE', 0*xE'];
% X0 = [ aph_b, 0, 0, 0, 0, 0];
X = fsolve( @(x) BarrierEndSystem( x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6) ), X0, options);

xOC = X(:,4)';
yOC = X(:,5)';
xC = X(:,2)';
yC = X(:,3)';
RC = hypot(xC, yC);
teC = atan2(yC, xC);


%% Outer base points D (bottom)
aphE = barrier_angles*pi/180;
teE = pi/2/p - aphE;
xE = Dend/2*cos(teE);
yE = Dend/2*sin(teE);

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

X0 = [ teE', xE', .8*yE', xE'*.9, yE'*.9, xE'*.1];
X = fsolve( @(x) BarrierEndSystem( x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6) ), X0, options);

xOD = X(:,4)';
yOD = X(:,5)';
xD = X(:,2)';
yD = X(:,3)';
RD = hypot(xD, yD);
teD = atan2(yD, xD);


%% Flux-barrier points
% We already have the potentials of the two flux-barrier sidelines
phiAprime = phi_fluid( rho_map(RAprime), xi_map(teAprime), rho0);
phiBprime = phi_fluid( rho_map(RBprime), xi_map(teBprime), rho0);

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

for bkk = 1:Nb
  dphiAC = (phiC(bkk) - phiAprime(bkk))./Nstep(bkk);
  dphiBD = (phiD(bkk) - phiBprime(bkk))./Nstep(bkk);
  % we create the matrix of potentials phi needed for points intersections
  PhiAC = phiAprime(bkk) + cumsum( repmat(dphiAC', Nstep(bkk) - 1, 1) );
  PhiBD = phiBprime(bkk) + cumsum( repmat(dphiBD', Nstep(bkk) - 1, 1) );
  PsiAC = repmat( psiA(bkk), Nstep(bkk)-1, 1);
  PsiBD = repmat( psiB(bkk), Nstep(bkk)-1, 1);
  
  X0 = [repmat(rho0*1.1, numel(PhiAC), 1), repmat(pi/4, numel(PhiAC), 1)];
  RhoXi_AC = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiAC, PhiAC, rho0 ), X0, options);
  RhoXi_BD = fsolve( @(x) PsiPhi( x(:,1),x(:,2), PsiBD, PhiBD, rho0 ), X0, options);
  
  R_AC = r_map(RhoXi_AC(:,1));
  te_AC = th_map(RhoXi_AC(:,2));
  R_BD = r_map(RhoXi_BD(:,1));
  te_BD = th_map(RhoXi_BD(:,2));
  
  if deb
    barrier(bkk).R_AC = R_AC;
    barrier(bkk).R_BD = R_BD;
    barrier(bkk).te_AC = te_AC;
    barrier(bkk).te_BD = te_BD;
  end
  
  % output of points
  barrier(bkk).Zeta = [...
    % top side
    xE(bkk) + j*yE(bkk)
    xOC(bkk) + j*yOC(bkk)
    xC(bkk) + j*yC(bkk)
    flipud( R_AC.*exp(j*te_AC) )
    RAprime(bkk).*exp(j*teAprime(bkk))
    % bottom side
    RBprime(bkk).*exp(j*teBprime(bkk))
    R_BD.*exp(j*te_BD)
    xD(bkk) + j*yD(bkk)
    xOD(bkk) + j*yOD(bkk)
    xE(bkk) + j*yE(bkk)
    ];
  
  barrier(bkk).X = real(barrier(bkk).Zeta);
  barrier(bkk).Y = imag(barrier(bkk).Zeta);
  
end


%% plot
if deb
  
  % draw the rotor
  figure
  hold on
  tt = linspace(0,pi/p,50);
  plot(R0*cos(tt), R0*sin(tt), 'k');
  plot(Dr/2*cos(tt), Dr/2*sin(tt), 'k');
  axis equal
  % plot the flux-barrier central point
  plot(RAprime.*exp(j*teAprime), 'rd')
  plot(RBprime.*exp(j*teBprime), 'bo')
  
  plot(xE, yE,'ko')
  
  plot(xOC, yOC,'go')
  plot(xC, yC,'ro')
  plot(xOD, yOD,'co')
  plot(xD, yD,'bo')
  
  %
  % plot(R_AC.*exp(j*te_AC),'r.-')
  % plot(R_BD.*exp(j*te_BD),'b.-')
  
  for bkk = 1:Nb
    plot(barrier(bkk).R_AC.*exp(j*barrier(bkk).te_AC),'r.-')
    plot(barrier(bkk).R_BD.*exp(j*barrier(bkk).te_BD),'b.-')
    
    plot(barrier(bkk).X, barrier(bkk).Y, '.-')
  end
  
end

end