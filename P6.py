#!/usr/bin/env python

# Predicting Particle Production in Power-Plant Plumes (P6) parametrizaton.
# Robin Stevens, Jan 2013

from pylab import *

# SVN revision number:
revstring = '$Rev:: 124  $'[6:-1]
#revnum = int(revstring)

# assumed values for unavailable parameters:
# (prefixed with ass for assumed, but also because you're an ass if
#  you are using them and you don't have to.)
# change the variable names and this comment to something nicer in the
# final product. 
assCS    = 0.01108 # background condensation sink [s-1]
assdswrf = 400.    # downward shortwave radiative flux [W/m2]
assvg    = 6.4     # boundary-layer wind speed [m/s]
assblh   = 500.    # boundary-layer height [m]
assd     = None    # distance from source [m]
assbgSO2 = 0.5     # background SO2 concentration [ppb]
assbgNOx = 1.0     # background NOx concentration [ppb]

minbgNOx = 5e-3 # minimum cutoff for bgNOx input [ppb]
minblh   = 50.  # minimum cutoff for boundary-layer height [m]

# Constants
eps       = 1e-10 # threshold for closure
dens      = 1770. # assuming (NH4)2SO4 [kg/m3]
N_Av      = 6.02214129e23 # Avagadro's constant [1/mol]
m_SO2     = 64.066e-3 # molar mass of SO4 [kg/mol]
m_H2SO4   = 98.08e-3 # molar mass of H2SO4 [kg/mol]
m_NH42SO4 = 132.14e-3 # molar mass of (NH4)2SO4 [kg/mol]
sig       = 1.4 # geometric standard deviation [unitless]



# Parameters for fitting to the mean concentration in the plume
# eg. bgNOx+p[0]*NOxemis*vg**p[1]*blh**p[2]*(time)**p[3]
so2fit = [1.70489332e+04,  -1.22925721e+00,  -1.89107567e-01,  -7.73243719e-01]
noxfit = [9.59537204e+04,  -1.23398130e+00,  -2.01833632e-01,  -7.90220955e-01]


def get_fox(NOxemis=None, dswrf=assdswrf, vg=assvg, blh=assblh, \
                d=assd,  bgNOx=assbgNOx, \
                p=[ -1.64966180e-10,   7.90402597e-01,   7.72321067e-01,   1.44390208e-08]):

    # p is the set of fitted parameters. This should never be given as
    # input, except by the other functions in this same module.

    locbgNOx = max(bgNOx,minbgNOx)
    time = d/vg

    effNOx = locbgNOx+p[3]*NOxemis \
        *vg**noxfit[1]*blh**noxfit[2]*(time)**noxfit[3]
    effNOx = max(effNOx,0.) # don't let this go less than zero

    # using Jim's parametrization in here:
    S0 = 1370. # Solar constant at top of atmosphere (W m-2)
    Tr = 0.76 # Transmittance through atmosphere
   
    p1p = [-0.014, 0.0027, 0.1713, -0.0466, -0.7893, -0.1739, 6.9414]
    p2p = [-1345, 4002, -471.8, 42.72]
   
    x = log10(effNOx)-0.195
    y = dswrf/(S0*Tr)
   
    p1 = p1p[0]*x**6+p1p[1]*x**5+p1p[2]*x**4+p1p[3]*x**3+p1p[4]*x**2+p1p[5]*x+p1p[6]
    p2 = (p2p[0]*y**3+p2p[1]*y**2+p2p[2]*y+p2p[3])*1.E4
    
    effOH = (0.82 * 10.**(p1*log10(p2)/6.8))
   
    fox = 1.-exp(p[0]*effOH**p[1]*time**p[2])

    assert (0.<=fox and fox<=1.)

    return fox

def get_nuc(SO2emis=None, NOxemis=None, CS=assCS, dswrf=assdswrf, \
                vg=assvg, blh=assblh, d=assd, bgSO2=assbgSO2, bgNOx=assbgNOx):

    p = [22387.211385683378, 1.92, 3.28, -1.24, -3.48, 436515.83224016655]
    locbgNOx = max(bgNOx,minbgNOx)
    time = d/vg

    if dswrf < 1e-3: # less than 1 mW/m2 - less than moonlight
        return False

    effSO2 = bgSO2+p[0]*SO2emis \
        *vg**so2fit[1]*blh**so2fit[2]*(time)**so2fit[3]
    if effSO2 < 1e-6: # less than 1e-3 ppt - extremely unlikely
        return False

    effNOx = bgNOx+p[5]*NOxemis \
        *vg**noxfit[1]*blh**noxfit[2]*(time)**noxfit[3]

    nucp = effSO2**p[1]*dswrf**p[2]*effNOx**p[3]*CS**p[4]+1e-10
    
    return nucp > 2.98841470581e+14

def get_mdiam(masspp):
    # for masspp in kg, returns number median diameter [um]

    mdiam = 1e6*(masspp/dens*6./pi)**(1./3.)*exp(-1.5*log(sig)**2.)

    assert mdiam >= 0.
    
    return mdiam

def get_masspp(SO2emis=None, NOxemis=None, CS=assCS, dswrf=assdswrf, \
                  vg=assvg, blh=assblh, d=assd, bgSO2=assbgSO2, bgNOx=assbgNOx):

    p = [ -1.29652905e-06,   6.92474330e-01,   2.92853444e-01,   2.13849343e+07,
           1.47496900e-27,   1.51723205e+00,   1.09357728e+00,  -6.17290992e-01,
           9.68490330e-01,   2.60502969e+06,   4.07112024e-23]

    time = d/vg

    locfox = get_fox(p=p[:4], NOxemis=NOxemis, dswrf=dswrf, vg=vg, \
                         blh=blh, d=d,  bgNOx=bgNOx)
    
    effSO2 = bgSO2 + p[9]*SO2emis*vg**so2fit[1]*blh**so2fit[2]*time**so2fit[3]
    masspp = p[4]*locfox**p[5]*effSO2**p[6]*CS**p[7]*time**p[8]+p[10]

    return masspp


def get_nnew(SO2emis=None, NOxemis=None, CS=assCS, dswrf=assdswrf, \
                  vg=assvg, blh=assblh, d=assd, bgSO2=assbgSO2, bgNOx=assbgNOx):

    p = [ -3.54855422e-15,  7.13304235e-01,    1.93747558e+00,   1.24321647e+06,
           6.93853928e+23,  9.94909098e-01,    2.49960504e-01,  -1.27968905e-01,
          -4.41706268e+00,  1.44126017e-01,    1.73637370e-01]


    time = d/vg

    locfox = get_fox(p=p[:4], NOxemis=NOxemis, dswrf=dswrf, vg=vg, \
                         blh=blh, d=d,  bgNOx=bgNOx)

    nnew = p[4] * locfox**p[5] * bgSO2**p[6] * SO2emis**p[7] \
        * exp( p[8] * CS**p[9] * time**p[10] ) + 1.

    assert nnew > 0.

    return nnew

def get_fnew(ifox=None,nuc=None,imasspp=None,innew=None):

    minmasspp = 2.*m_H2SO4/N_Av # mass of two SO4 molecules

    ffox = ifox
    fmasspp = imasspp
    fnnew = innew

    if nuc:
    
        ifnew = (fmasspp*fnnew)/ffox*m_SO2/m_H2SO4
        ffnew = ifnew

        if ffnew > 1.:
            # where cfnew > 1., reduce it to one, and split the difference
            # among masspp and nnew for closure.
            fmasspp /= ffnew**(0.5)
            fnnew   /= ffnew**(0.5)
            # I don't think this condition is possible anymore,
            # but just in case...
            if fmasspp < minmasspp: 
                fnnew *= fmasspp/minmasspp
                fmasspp = minmasspp
            ffnew = 1.
        if abs((fmasspp*fnnew)/(ffox*ffnew)*m_SO2/m_H2SO4 -1.)>eps:
            print('closure not achieved in get_fnew')
            print('fnew',ffnew)
            print('current ratio',(fmasspp*fnnew)/(ffox*ffnew)*m_SO2/m_H2SO4)

    else:
        ffnew = 0.

    assert (0.<=ffox and ffox<=1.)
    assert fmasspp > minmasspp
    assert fnnew > 0.
    assert (0.<=ffnew and ffnew<=1.)

    return ffox,fmasspp,fnnew,ffnew

def get_em_all(SO2emis=None, NOxemis=None, CS=assCS, \
                   dswrf=assdswrf, vg=assvg, blh=assblh, \
                   d=assd, bgSO2=assbgSO2, bgNOx=assbgNOx):

    blh = max(blh, minblh)

    fox = get_fox(NOxemis=NOxemis, dswrf=dswrf, vg=vg, \
                       blh=blh, d=d, bgNOx=bgNOx)


    nuc    =   get_nuc(SO2emis=SO2emis, NOxemis=NOxemis, CS=CS, \
                           dswrf=dswrf, vg=vg, blh=blh, \
                           d=d, bgSO2=bgSO2, bgNOx=bgNOx)


    masspp = get_masspp(SO2emis=SO2emis, NOxemis=NOxemis, CS=CS, \
                           dswrf=dswrf, vg=vg, blh=blh, \
                           d=d, bgSO2=bgSO2, bgNOx=bgNOx)
    

    nnew  =  get_nnew(SO2emis=SO2emis, NOxemis=NOxemis, CS=CS, \
                           dswrf=dswrf, vg=vg, blh=blh, \
                           d=d, bgSO2=bgSO2, bgNOx=bgNOx)

    fox,masspp,nnew,fnew  =  get_fnew(fox,nuc,masspp,nnew)


    return(fox,nuc,masspp,nnew,fnew)

def splitemis(id=assd,iSO2emis=None,iNOxemis=None,iCS=assCS,iDSWRF=assdswrf, \
                  ivg=assvg,iBLH=assblh,ibgSO2=assbgSO2,ibgNOx=assbgNOx):

    # This subroutine calls P6 using an assumed distribution of low-,
    # medium-, and high-emitting sources

    #-----ADJUSTABLE PARAMETERS-------------------------------------------

    # default values for SO2 emissions:
    lSO2emis=0.0606 # low  SO2 emissions [kg/s]
    mSO2emis=0.202  # med  SO2 emissions [kg/s]
    hSO2emis=1.00   # high SO2 emissions [kg/s]
  
    # default values for NOx emissions:
    dlNOxemis=0.0300# low  NOx emissions [kg N/s]
    dmNOxemis=0.0840# med  NOx emissions [kg N/s]
    dhNOxemis=0.290 # high NOx emissions [kg N/s]

    #-----CODE--------------------------------------------------------------
      
    if iSO2emis == None or iNOxemis == None:
        lNOxemis = dlNOxemis
        mNOxemis = dmNOxemis
        hNOxemis = dhNOxemis
    else:
        lNOxemis = iNOxemis/iSO2emis * lSO2emis
        mNOxemis = iNOxemis/iSO2emis * mSO2emis
        hNOxemis = iNOxemis/iSO2emis * hSO2emis
    

    lfox,lnuc,lmasspp,lnnew,lfnew = \
        get_em_all(SO2emis=lSO2emis, NOxemis=lNOxemis, CS=iCS, \
                       dswrf=iDSWRF, vg=ivg, blh=iBLH, \
                       d=id, bgSO2=ibgSO2, bgNOx=ibgNOx)

    mfox,mnuc,mmasspp,mnnew,mfnew = \
        get_em_all(SO2emis=mSO2emis, NOxemis=mNOxemis, CS=iCS, \
                       dswrf=iDSWRF, vg=ivg, blh=iBLH, \
                       d=id, bgSO2=ibgSO2, bgNOx=ibgNOx)

    hfox,hnuc,hmasspp,hnnew,hfnew = \
        get_em_all(SO2emis=hSO2emis, NOxemis=hNOxemis, CS=iCS, \
                       dswrf=iDSWRF, vg=ivg, blh=iBLH, \
                       d=id, bgSO2=ibgSO2, bgNOx=ibgNOx)    

    # overall fox is emissions-weighted average of lfox, mfox, and hfox
    ofox = (lfox*lSO2emis + mfox*mSO2emis + hfox*hSO2emis)/ \
        (lSO2emis+mSO2emis+hSO2emis)

    if (lnuc==0. and mnuc==0. and hnuc==0.) :
        # no nucleation
        omasspp = 0.
        omdiam = 0.
        onnew = 0.
        ofnew = 0.
    else:
        # overall mass per particle is average weighted by number of 
        # particles, disregarding no nucleation cases
        omasspp = (lmasspp*lSO2emis*lnnew*lnuc + \
                       mmasspp*mSO2emis*mnnew*mnuc + \
                       hmasspp*hSO2emis*hnnew*hnuc)/ \
                  (lSO2emis*lnnew*lnuc+mSO2emis*mnnew*mnuc+hSO2emis*hnnew*hnuc)
        # overall nnew is total new particles over total emissions
        onnew = (lnnew*lSO2emis + mnnew*mSO2emis + hnnew*hSO2emis)/ \
            (lSO2emis+mSO2emis+hSO2emis)

        # overall fnew is from closure
        # note that fox*SO2emis weighted average would give the same answer
        tmpmasspp = omasspp
        tmpnnew = onnew
        ofox,omasspp,onnew,ofnew = get_fnew(ofox,1,tmpmasspp,tmpnnew)
        
        # overall mdiam just calculated from omasspp
        omdiam = get_mdiam(omasspp)

    omdiam = get_mdiam(omasspp)

    return(ofox,omasspp,omdiam,onnew,ofnew)

