# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 21:49:54 2018

@author: Giacomo
"""

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

# I define a dummy class for Matlab structure like objects
class structtype():
    pass


def psi_fluid(rho,xi,rho0):
    return (rho**2 - rho0**2)/rho*np.sin(xi);
#    return (np.power(rho,2) - rho0**2)/rho*np.sin(xi);

def phi_fluid(rho,xi,rho0):
    return (np.power(rho,2) + rho0**2)/rho*np.cos(xi);

def xi_fluid(psi,rho,rho0):
    return np.arcsin(psi*rho/(np.power(rho,2) - rho0**2));

def rho_fluid(psi,xi,rho0):
    return ( psi + np.sqrt(np.power(psi,2) + 4*np.power(np.sin(xi),2)*rho0**2) )/(2*np.sin(xi));

def r_map(rho):
    return np.power(rho, 1/p);

def th_map(xi):
    return xi/p;

def rho_map(r):
    return np.power(r, p);

def xi_map(th):
    return th*p;
    
def vr(r,th,R0):
    return p*(np.power(r,(p-1)) - R0**(2*p)/np.power(r,(p+1)))*np.cos(p*th);

def vt(r,th,R0):
    return -p*( np.power(r,(p-1)) + R0**(2*p)/np.power(r,(p+1)) )*np.sin(p*th);

def vx(vr_v,vth_v,th):
    return vr_v*np.cos(th) - vth_v*np.sin(th);

def vy(vr_v,vth_v,th):
    return vr_v*np.sin(th) + vth_v*np.cos(th);


def CentralPt_Eq(th, *args):
    psiCentralPt, rho0, mCentral, qCentral = args;
    return np.multiply( r_map( rho_fluid(psiCentralPt, xi_map(th), rho0) ),
                 ( np.sin(th) - mCentral*np.cos(th) ) ) - qCentral;


def BarrierEndSystem(X, *args):
#    th,xd,yd,xo,yo,R = X;
    th = X[0:Nb];
    xd = X[1*Nb:2*Nb];
    yd = X[2*Nb:3*Nb];
    xo = X[3*Nb:4*Nb];
    yo = X[4*Nb:5*Nb];
    R  = X[5*Nb:6*Nb];
    
    psiA, rho0, xE, yE = args;
    R0 = r_map(rho0);
    firstEq = xd - np.multiply( r_map(rho_fluid(psiA, p*th, rho0)), np.cos(th) );
    seconEq = yd - np.multiply( r_map(rho_fluid(psiA, p*th, rho0)), np.sin(th) );
    thirdEq = (xd - xo)**2 + (yd - yo)**2 - R**2;
#    thirdEq = (xE - xo)**2 + (yE - yo)**2 - R**2;
    fourtEq = (xE - xo)**2 + (yE - yo)**2 - R**2;
    fifthEq = np.multiply( (xo - xd), vx( vr( r_map(rho_fluid(psiA, p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiA, p*th, rho0)) ,th,R0 ), th) ) + np.multiply( (yo - yd), vy( vr( r_map(rho_fluid(psiA, p*th, rho0)),th,R0 ), vt( r_map(rho_fluid(psiA, p*th, rho0)) ,th,R0 ), th) );
    sixthEq = np.multiply(xo - xE,yE) - np.multiply(yo - yE,xE);
    return np.concatenate([firstEq,
                           seconEq,
                           thirdEq,
                           fourtEq,
                           fifthEq,
                           sixthEq])

def PsiPhi(X, *args):
    psi, phi, rho0, N = args;
    rho = X[0:N];
    xi  = X[N:2*N];
    return np.concatenate( [psi - psi_fluid(rho, xi, rho0),
            phi - phi_fluid(rho, xi, rho0)] )


# MAIN function definition
def calc_fluid_barrier(r, deb):
    "CALC_FLUID_BARRIER computes the flux-barrier points along the streamline function."
    
    Dr = r.De; # [m], rotor outer diameter
    ScalingFactor = 1/( 10**(round(np.log10(Dr))) );
    # ScalingFactor = 1;
    Dr = Dr*ScalingFactor;
    
    pi = np.pi;
    global p, Nb # I have been lazy here...
    p = r.p; # number of pole pairs
    Nb = r.Nb; # number of flux-barriers
    tb = r.tb*ScalingFactor; # flux-barrier widths
    wc = r.wc*ScalingFactor; # flux-carrier widths
    Nstep = r.Nstep; # number of steps to draw the flux-barrier side
    
    wrib_t = r.wrib_t*ScalingFactor; # [m], tangential iron rib width

    
    if hasattr(r,'barrier_angles_el'):
        barrier_angles_el = r.barrier_angles_el; # [deg], electrical flux-barrier angles
        AutoBarrierEndCalc = 0;
    else:
        barrier_angles_el =  np.zeros([1,Nb]);
        AutoBarrierEndCalc = 1;
    
    if hasattr(r,'wm'):
        wm = r.wm*ScalingFactor;
    else:
        wm = 0;
      
    if hasattr(r,'wrib'):
        wrib = r.wrib*ScalingFactor + wm; # [m], radial iron rib widths
    else:
        wrib = np.zeros([1,Nb]) + wm;
    
    Dend = Dr - 2*wrib_t; # [m], flux-barrier end diameter
    Dsh = Dend - 2*( np.sum(tb) + np.sum(wc) ); # [m], shaft diameter
    R0 = Dsh/2; # [m], shaft radius
    barrier_angles = barrier_angles_el/p; # [deg], flux-barrier angles
    if hasattr(r,'barrier_end'):
        barrier_end = r.barrier_end;
    else:
        barrier_end = '';

    ## Precomputations
    rho0 = rho_map(R0);

    ## Central base points 
    RAprime = Dend/2 - np.concatenate( (np.array([0]), np.cumsum( tb[0:-1]) ) ) - np.cumsum( wc[0:-1] ); # top
    RBprime = RAprime - tb; # bottom
    te_qAxis = pi/(2*p); # q-axis angle in rotor reference frame
    
    # get A' and B' considering rib and magnet widths
    mCentral = np.tan(te_qAxis); # slope
    qCentral = np.tile( -wrib/2/np.cos(te_qAxis), 2); # intercept
    
    psiCentralPtA = psi_fluid(rho_map(RAprime), xi_map(te_qAxis), rho0);
    psiCentralPtB = psi_fluid(rho_map(RBprime), xi_map(te_qAxis), rho0);
    psiCentralPt = np.array( np.concatenate( (psiCentralPtA, psiCentralPtB) ) );
    psiA = psiCentralPtA;
    psiB = psiCentralPtB;
    

    FunctionTolerance = 10*np.spacing(1);
    StepTolerance = 1e4*np.spacing(1);

    X0 = np.repeat(te_qAxis, 2*Nb);
    data = (psiCentralPt,rho0,mCentral,qCentral);
    # test function 
#    print( CentralPt_Eq(X0, *data ) )

    teAB = sp.optimize.fsolve(CentralPt_Eq, X0, args=data, xtol=StepTolerance, epsfcn=FunctionTolerance);
    teA = teAB[0:Nb];
    teB = teAB[Nb:];
    RA = r_map( rho_fluid(psiA, xi_map(teA), rho0) );
    RB = r_map( rho_fluid(psiB, xi_map(teB), rho0) );
    
    # central base points
    zA = RA*np.exp(1j*teA);
    zB = RB*np.exp(1j*teB);
    xA = zA.real;
    yA = zA.imag;
    xB = zB.real;
    yB = zB.imag;
    
    # magnet central base point radius computation
    RAsecond = RA*np.cos(te_qAxis - teA);
    RBsecond = RB*np.cos(te_qAxis - teB);  
    
    Rmag = (RAprime + RAsecond + RBprime + RBsecond)/4;
    
    # 1st test --> OK!
#    print(RA,teA,RB,teB)
#    print(xA,yA,xB,yB)

    # Outer base points C,D preparation
    RCprime = Dend/2;
    teCprime = th_map( xi_fluid(psiA, rho_map(RCprime), rho0) );
    xCprime = Dend/2*np.cos(teCprime);
    yCprime = Dend/2*np.sin(teCprime);
    
    RDprime = Dend/2;
    teDprime = th_map( xi_fluid(psiB, rho_map(RDprime), rho0) );
    xDprime = Dend/2*np.cos(teDprime);
    yDprime = Dend/2*np.sin(teDprime);

    if AutoBarrierEndCalc:
        teE = (teCprime + teDprime)/2;
        aphE = pi/2/p - teE;
        barrier_angles = 180/np.pi*aphE;
        barrier_angles_el = p*barrier_angles;
    else:
        aphE = barrier_angles*pi/180;
        teE = pi/2/p - aphE;

    xE = Dend/2*np.cos(teE);
    yE = Dend/2*np.sin(teE);  
    
    # 2nd test --> OK!
#    print(xE,yE)
    










    ## Outer base points C (top)
    if barrier_end == 'rect':
        RC = RCprime;
        teC = teCprime;
        xC = xCprime;
        yC = yCprime;
        xOC = xC;
        yOC = yC;
    
    else:      
        # 1st try
#        X0 = [ 1.5*teE, 0.9*xE, 0.9*yE, 0.8*xE, 0.8*yE, 0.25*xE];
        # best try
        xC0 = ( xE + xCprime + 0.1*xA )/(2 + 0.1);
        yC0 = ( yE + yCprime )/2;
        thC0 = np.arctan(yC0/xC0);
        xOC0 = ( xE + xC0 + 0 )/3;
        yOC0 = ( yE + yC0 + 0 )/3;
        RCOCE0 = np.sqrt( (xOC0 - xE)**2 + (yOC0 - yE)**2 );
        
        X0 = [ thC0, xC0, yC0, xOC0, yOC0, RCOCE0];
        X0 = np.reshape(X0, Nb*6);
      
        data = (psiA, rho0, xE, yE);
        X = sp.optimize.fsolve( BarrierEndSystem, X0, args=data);
      
        xOC = X[3*Nb:4*Nb];
        yOC = X[4*Nb:5*Nb];
        xC = X[1*Nb:2*Nb];
        yC = X[2*Nb:3*Nb];
        RC = np.sqrt(xC**2 + yC**2);
        teC = np.arctan2(yC, xC);
        
    # 3rd test --> OK!
#    print(xOC)
#    print(yOC)
#    print(xC)
#    print(yC)
#    print(RC)
#    print(teC)


    ## Outer base points D (bottom)
    if barrier_end == 'rect':
        RD = RDprime;
        teD = teDprime;
        xD = xDprime;
        yD = yDprime;
        xOD = xD;
        yOD = yD;
    
    else:      
        # 1st try
#        X0 = [ 0.8*teE, 0.8*xE, 0.8*yE, 0.9*xE, 0.9*yE, 0.2*xE];
        # best try
        xD0 = ( xE + xDprime )/2;
        yD0 = ( yE + yDprime )/2;
        thD0 = np.arctan(yD0/xD0);
        xOD0 = ( xE + xD0 + xC )/3;
        yOD0 = ( yE + yD0 + yC )/3;
        RDODE0 = np.sqrt( (xOD0 - xE)**2 + (yOD0 - yE)**2 );
        
        X0 = [ thD0, xD0, yD0, xOD0, yOD0, RDODE0];
        X0 = np.reshape(X0, Nb*6);
      
        data = (psiB, rho0, xE, yE);
        X = sp.optimize.fsolve( BarrierEndSystem, X0, args=data);
      
        xOD = X[3*Nb:4*Nb];
        yOD = X[4*Nb:5*Nb];
        xD = X[1*Nb:2*Nb];
        yD = X[2*Nb:3*Nb];
        RD = np.sqrt(xD**2 + yD**2);
        teD = np.arctan2(yD, xD);
        
    # 4th test --> OK!
#    print(xOD)
#    print(yOD)
#    print(xD)
#    print(yD)
#    print(RD)
#    print(teD)


    ## Flux-barrier points
    # We already have the potentials of the two flux-barrier sidelines
    phiA = phi_fluid( rho_map(RA), xi_map(teA), rho0);
    phiB = phi_fluid( rho_map(RB), xi_map(teB), rho0);
    
    phiC = phi_fluid( rho_map(RC), xi_map(teC), rho0);
    phiD = phi_fluid( rho_map(RD), xi_map(teD), rho0);
    
    barrier = structtype();

    
    
    
    
    
    
    
    
    
    
    
    XX = [];
    YY = [];
    Rm = [];
    
    for bkk in range(0,Nb):
        dphiAC = np.divide(phiC[bkk] - phiA[bkk], Nstep[bkk]);
        dphiBD = np.divide(phiD[bkk] - phiB[bkk], Nstep[bkk]);
        # we create the matrix of potentials phi needed for points intersections
        PhiAC = phiA[bkk] + np.cumsum( np.tile(dphiAC, Nstep[bkk] - 1) );
        PhiBD = phiB[bkk] + np.cumsum( np.tile(dphiBD, Nstep[bkk] - 1) );
        PsiAC = np.tile( psiA[bkk], Nstep[bkk]-1);
        PsiBD = np.tile( psiB[bkk], Nstep[bkk]-1);

        X0 = np.concatenate( [np.linspace(rho0, Dend/2, np.size(PhiAC)), np.linspace(pi/4, xi_map(teE[bkk]), np.size(PhiAC))] );

        data = (PsiAC, PhiAC, rho0, Nstep[bkk]-1);
        RhoXi_AC = sp.optimize.fsolve( PsiPhi, X0, args=data );

        data = (PsiBD, PhiBD, rho0, Nstep[bkk]-1);
        RhoXi_BD = sp.optimize.fsolve( PsiPhi, X0, args=data );
    
        R_AC = r_map( RhoXi_AC[0:Nstep[bkk]-1] );
        te_AC = th_map( RhoXi_AC[Nstep[bkk]-1:] );
        R_BD = r_map( RhoXi_BD[0:Nstep[bkk]-1] );
        te_BD = th_map( RhoXi_BD[Nstep[bkk]-1:] );
        
        # 5th test --> OK!
#        print(R_AC, te_AC)
#        print(R_BD, te_BD)
          
        Zeta = np.concatenate( [[
                # top side
                xE[bkk] + 1j*yE[bkk],
                xOC[bkk] + 1j*yOC[bkk],
                xC[bkk] + 1j*yC[bkk] ],
                np.flipud( np.multiply(R_AC, np.exp(1j*te_AC)) ),
                [xA[bkk] + 1j*yA[bkk],
                 # bottom  side
                 xB[bkk] + 1j*yB[bkk]],
                np.multiply( R_BD, np.exp(1j*te_BD) ),                
                [xD[bkk] + 1j*yD[bkk],
                 xOD[bkk] + 1j*yOD[bkk],
                 xE[bkk] + 1j*yE[bkk]]           
        ] )/ScalingFactor;

        X = Zeta.real;
        Y = Zeta.imag;
        
        XX.append([X]);
        YY.append([Y]);
        
        # magnet central base point
        Rm.append(Rmag[bkk]/ScalingFactor);

    barrier.X = XX;
    barrier.Y = YY;
    barrier.Rm = Rm;



    if deb == 1:
        Apt = np.multiply(RA/ScalingFactor, np.exp(1j*teA) );
        Bpt = np.multiply(RB/ScalingFactor, np.exp(1j*teB) );
        Cpt = np.multiply(RC/ScalingFactor, np.exp(1j*teC) );
        Dpt = np.multiply(RD/ScalingFactor, np.exp(1j*teD) );
        plt.plot( [Apt.real,Cpt.real], [Apt.imag,Cpt.imag], "ro" )
        plt.plot( [Bpt.real,Dpt.real], [Bpt.imag,Dpt.imag], "bo" )
        plt.plot( xE/ScalingFactor, yE/ScalingFactor, "ko" )
        
        for bkk in range(0,Nb):
            plt.plot(np.squeeze(barrier.X[bkk]), np.squeeze(barrier.Y[bkk]))
            
        plt.axis('equal')
        plt.show()
    
        
    
    
    
    return barrier