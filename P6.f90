! **************************************************
! *  P6                                            *
! **************************************************

! Written by Robin Stevens, June 2013

! Purpose: The Predicting Particles Produced in Power-Plant Plumes
! parameterization calculates the characteristics of aerosol formed
! in sulphur-rich plumes based on the emissions, meteorological
! conditions, background aersol condensation sink and background SO2
! and NOx concentrations.

! The input parameters and their tested ranges are described in the
! P6 "INPUTS" section. The subroutine should give sensible results
! outside of the tested range, but these values have not been tested.

! The subroutine is described in the accompanying paper:
! R. G. Stevens and J. R. Pierce,
! A parameterization of sub-grid particle formation in sulphur-rich
! plumes for global and regional-scale models.,
! Atmos. Chem. Phys. Discuss., 13, 19583-19623, 2013
! www.atmos-chem-phys-discuss.net/13/19583/2013/
! doi:10.5194/acpd-13-19583-2013

!-----USAGE-----------------------------------------------------------

! For a single source, with all inputs available, call P6 with each
! input and output variable specified.

! Examples:
! call P6(id=d,iSO2emis=SO2emis,iNOxemis=NOxemis,iCS=CS, &
!     iDSWRF=DSWRF,ivg=vg,iBLH=BLH,ibgSO2=bgSO2,ibgNOx=bgNOx, &
!     ofox=fox,omasspp=masspp,omdiam=mdiam,onnew=nnew,ofnew=fnew)
! call P6(id,iSO2emis,iNOxemis,iCS,iDSWRF,ivg,iBLH,ibgSO2,ibgNOx, &
!     ofox,omasspp,omdiam,onnew,ofnew)

! For a single source, with some inputs available, call P6 with only
! the available inputs and the output variables specified. A
! pre-defined value will be assumed for the missing inputs, listed in
! the adjustable parameters for P6. If NOxemis is missing, a typical
! NOx:SO2 emissions ratio will be assumed. Note that id and iSO2emis 
! are non-optional. Note that you must specify each argument by name
! if any are missing.

! Examples:
! call P6(id=d,iSO2emis=SO2emis,iNOxemis=NOxemis,iCS=CS, &
!     ibgSO2=bgSO2,ibgNOx=bgNOx, &
!     ofox=fox,omasspp=masspp,omdiam=mdiam,onnew=nnew,ofnew=fnew)
! call P6(id=d,iSO2emis=SO2emis, &
!     ofox=fox,omasspp=masspp,omdiam=mdiam,onnew=nnew,ofnew=fnew)

! If only grid-based emissions are available, or if SO2 emissions are
! not available, call P6_splitemis to split the emissions amongst a
! set of sources based on 2010 US coal-fired power plant emissions.
! If both iSO2emis and iNOxemis are given, the distribution of SO2
! emissions will be assumed, and the NOx:SO2 ratio will be the given
! iNOxemis:iSO2emis ratio. If either iSO2emis or iNOxemis is missing,
! both the SO2emissions and the NOxemissions will be assumed.
! The only input that is non-optional is id.

! Examples:
! call P6_splitemis(id=d,iSO2emis=SO2emis,iNOxemis=NOxemis,iCS=CS, &
!     iDSWRF=DSWRF,ivg=vg,iBLH=BLH,ibgSO2=bgSO2,ibgNOx=bgNOx, &
!     ofox=fox,omasspp=masspp,omdiam=mdiam,onnew=nnew,ofnew=fnew)
! call P6_splitemis(id,iSO2emis,iNOxemis,iCS,iDSWRF,ivg,iBLH,ibgSO2,ibgNOx, &
!     ofox,omasspp,omdiam,onnew,ofnew)
! call P6_splitemis(id=d, &
!     ofox=fox,omasspp=masspp,omdiam=mdiam,onnew=nnew,ofnew=fnew)

module P6mod

implicit none


!-----MODULE PARAMETERS-----------------------------------------------

double precision, parameter:: eps = 1d-10 ! a very small number

double precision so2fit(3) ! fitted parameters for effective in-plume SO2
double precision noxfit(3) ! fitted parameters for effective in-plume NOx
data so2fit / -1.22925721d+00,  -1.89107567d-01,  -7.73243719d-01 /
data noxfit / -1.23398130d+00,  -2.01833632d-01,  -7.90220955d-01 /

contains

  subroutine P6(id,iSO2emis,iNOxemis,iCS,iDSWRF,ivg,iBLH,ibgSO2,ibgNOx, &
     ofox,omasspp,omdiam,onnew,ofnew)
    ! Returns the characteristics of sub-grid sulphate formed within a
    ! sulphur-rich plume based on the emissions, meteorological conditions,
    ! background aersol condensation sink and background SO2 and NOx.

    !-----INPUTS------------------------------------------------------------

    double precision, intent(in):: id
    ! distance from source [m]
    ! Tested range: 5000 - 100000 m
    double precision, intent(in):: iSO2emis
    ! SO2 emissions [kg/s]
    ! Tested range: 0.001 - 10 kg/s
    double precision, intent(in), optional:: iNOxemis
    ! NOx emissions [kg N/s]
    ! Tested range: 0.001 - 2 kg/s       pre-defined NOx:SO2 ratio: 0.419
    double precision, intent(in), optional:: iCS
    ! background condensation sink [s-1]
    ! Tested range: 8.94e-5 - 1.46e-2 s-1 pre-defined value: 1.108e-2 s-1
    double precision, intent(in), optional:: iDSWRF
    ! downward shortwave radiative flux [W/m2]
    ! Tested range: 100 - 960 W/m2       pre-defined value: 400 W/m2
    double precision, intent(in), optional:: ivg
    ! mean boundary-layer wind speed [m/s]
    ! Tested range: 0.178 - 26.1 m/s     pre-defined value: 6.4 m/s
    double precision, intent(in), optional:: iBLH
    ! boundary layer height [m]
    ! Tested range: 53 - 2792 m          pre-defined value: 500 m
    double precision, intent(in), optional:: ibgSO2
    ! background SO2 concentration [ppb]
    ! Tested range: 1.27e-6 - 16.6 ppb   pre-defined value: 0.5 ppb
    double precision, intent(in), optional:: ibgNOx
    ! background NOx concentration [ppb]
    ! Tested range: 2.84e-4 - 7.93 ppb   pre-defined value: 1.0 ppb

    !-----OUTPUTS-----------------------------------------------------------

    double precision, intent(out):: ofox
    ! fraction of emitted SO2 oxidized [unitless]
    double precision, intent(out):: omasspp
    ! mass per particle of newly-formed sulphate particles [kg]
    double precision, intent(out):: omdiam
    ! median diameter of newly-formed sulphate particles [nm]
    double precision, intent(out):: onnew
    ! number of new particles per kg SO2 emitted [kg-1 SO2]
    double precision, intent(out):: ofnew
    ! fraction of H2SO4 formed that forms new particles [unitless]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    ! the inputs, again:
    double precision SO2emis, NOxemis, DSWRF, CS, vg, BLH, d, bgSO2, bgNOx
    double precision tmpmasspp
    double precision tmpnnew

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    ! default values for unavailable inputs:
    double precision, parameter:: defCS    = 0.01108d0 ! condensation sink [s-1]
    double precision, parameter:: defDSWRF = 400.d0    ! downward shortwave radiative flux [W/m2]
    double precision, parameter:: defvg    = 6.4d0     ! boundary layer wind speed [m/s]
    double precision, parameter:: defBLH   = 500.d0    ! boundary layer height [m]
    double precision, parameter:: defbgSO2 = 0.5d0     ! background SO2 concentration [ppb]
    double precision, parameter:: defbgNOx = 1.0d0     ! background NOx concentration [ppb]
    
    double precision, parameter:: noxso2rat = 0.419d0  ! assumed NOx:SO2 ratio, if NOxemis not given
    double precision, parameter:: minbgNOx = 5d-3      ! minimum cutoff for bgNOx input [ppb]


    !-----CODE--------------------------------------------------------------

    d = id

    SO2emis = iSO2emis

    if (present(iNOxemis)) then
       NOxemis = iNOxemis
    else
       NOxemis = SO2emis * noxso2rat
    endif

    if (present(iCS)) then
       CS = iCS
    else
       CS = defCS
    endif

    if (present(iDSWRF)) then
       DSWRF = iDSWRF
    else
       DSWRF = defDSWRF
    endif

    if (present(ivg)) then
       vg = ivg
    else
       vg = defvg
    endif

    if (present(iBLH)) then
       BLH = iBLH
    else
       BLH = defBLH
    endif

    if (present(ibgSO2)) then
       bgSO2 = ibgSO2
    else
       bgSO2 = defbgSO2
    endif

    if (present(ibgNOx)) then
       bgNOx = max(ibgNOx,minbgNOx)
    else
       bgNOx = defbgNOx
    endif

    ofox = getfox(d,NOxemis,DSWRF,vg,BLH,bgNOx)

    if (getnuc(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)) then
       ! not-insignificant new-particle formation

       omasspp = getmasspp(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)
       onnew = getnnew(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)

       tmpmasspp = omasspp
       tmpnnew = onnew
       call getfnew(ofox,tmpmasspp,tmpnnew,omasspp,onnew,ofnew)

       omdiam = getmdiam(omasspp)

    else
       ! insignificant new-particle formation     
       omasspp = 0.
       omdiam = 0.
       onnew = 0.
       ofnew = 0.
       
    endif


  end subroutine P6


  subroutine P6_splitemis(id,iSO2emis,iNOxemis,iCS,iDSWRF, &
       ivg,iBLH,ibgSO2,ibgNOx, &
       ofox,omasspp,omdiam,onnew,ofnew)
    ! This subroutine calls P6 using an assumed distribution of low-,
    ! medium-, and high-emitting sources

    !-----INPUTS------------------------------------------------------------

    double precision, intent(in):: id
    ! distance from source [m]
    double precision, intent(in), optional:: iSO2emis
    ! SO2 emissions [kg/s]
    double precision, intent(in), optional:: iNOxemis
    ! NOx emissions [kg N/s]
    double precision, intent(in), optional:: iCS 
    ! background condensation sink [s-1]
    double precision, intent(in), optional:: iDSWRF
    ! downward shortwave radiative flux [W/m2]
    double precision, intent(in), optional:: ivg
    ! mean boundary-layer wind speed [m/s]
    double precision, intent(in), optional:: iBLH
    ! boundary layer height [m]
    double precision, intent(in), optional:: ibgSO2
    ! background SO2 concentration [ppb]
    double precision, intent(in), optional:: ibgNOx
    ! background NOx concentration [ppb]
    
    !-----OUTPUTS-----------------------------------------------------------
    
    double precision, intent(out):: ofox
    ! fraction of emitted SO2 oxidized [unitless]
    double precision, intent(out):: omasspp 
    ! mass per particle of newly-formed sulphate particles [kg]
    double precision, intent(out):: omdiam
    ! median diameter of newly-formed sulphate particles [nm]
    double precision, intent(out):: onnew 
    ! number of new particles per kg SO2 emitted [kg-1 SO2]
    double precision, intent(out):: ofnew 
    ! fraction of H2SO4 formed that forms new particles [unitless]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    double precision lnuc, mnuc, hnuc ! is there nucleation in each case?

    double precision nsources ! the number (each) of low, med, and high emitters

    double precision tmpmasspp
    double precision tmpnnew

    double precision lNOxemis ! low  NOx emissions [kg/s]
    double precision mNOxemis ! med  NOx emissions [kg/s]
    double precision hNOxemis ! high NOx emissions [kg/s]

    ! outputs from high, low, and medium emitters:
    double precision lfox, lmasspp, lmdiam, lnnew, lfnew
    double precision mfox, mmasspp, mmdiam, mnnew, mfnew
    double precision hfox, hmasspp, hmdiam, hnnew, hfnew

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    double precision, parameter:: lSO2emis=0.0606d0 ! low  SO2 emissions [kg/s]
    double precision, parameter:: mSO2emis=0.202d0  ! med  SO2 emissions [kg/s]
    double precision, parameter:: hSO2emis=1.00d0   ! high SO2 emissions [kg/s]
  
    ! default values for NOx emissions:
    double precision, parameter:: dlNOxemis=0.0300d0! low  NOx emissions [kg/s]
    double precision, parameter:: dmNOxemis=0.0840d0! med  NOx emissions [kg/s]
    double precision, parameter:: dhNOxemis=0.290d0 ! high NOx emissions [kg/s]

    !-----CODE--------------------------------------------------------------
      
    if (present(iNOxemis) .and. present(iSO2emis)) then
       lNOxemis = iNOxemis/iSO2emis * lSO2emis
       mNOxemis = iNOxemis/iSO2emis * mSO2emis
       hNOxemis = iNOxemis/iSO2emis * hSO2emis
    else
       lNOxemis = dlNOxemis
       mNOxemis = dmNOxemis
       hNOxemis = dhNOxemis
    endif

    call P6(id=id,iSO2emis=lSO2emis,iNOxemis=lNOxemis,iDSWRF=iDSWRF, &
         iCS=iCS,ivg=ivg,iBLH=iBLH,ibgSO2=ibgSO2,ibgNOx=ibgNOx, &
         ofox=lfox,omasspp=lmasspp,omdiam=lmdiam,onnew=lnnew,ofnew=lfnew)

    ! nnew = 0. and fnew = 0. implies no nucleation
    if (lnnew .lt. eps .and. lfnew .lt. eps) then
       lnuc = 0.d0
    else
       lnuc = 1.d0
    endif


    call P6(id=id,iSO2emis=mSO2emis,iNOxemis=mNOxemis,iDSWRF=iDSWRF, &
         iCS=iCS,ivg=ivg,iBLH=iBLH,ibgSO2=ibgSO2,ibgNOx=ibgNOx, &
         ofox=mfox,omasspp=mmasspp,omdiam=mmdiam,onnew=mnnew,ofnew=mfnew)

    ! nnew = 0. and fnew = 0. implies no nucleation
    if (mnnew .lt. eps .and. mfnew .lt. eps) then
       mnuc = 0.d0
    else
       mnuc = 1.d0
    endif

    call P6(id=id,iSO2emis=hSO2emis,iNOxemis=hNOxemis,iDSWRF=iDSWRF, &
         iCS=iCS,ivg=ivg,iBLH=iBLH,ibgSO2=ibgSO2,ibgNOx=ibgNOx, &
         ofox=hfox,omasspp=hmasspp,omdiam=hmdiam,onnew=hnnew,ofnew=hfnew)

    ! nnew = 0. and fnew = 0. implies no nucleation
    if (hnnew .lt. eps .and. hfnew .lt. eps) then
       hnuc = 0.d0
    else
       hnuc = 1.d0
    endif

    ! overall fox is emissions-weighted average of lfox, mfox, and hfox
    ofox = (lfox*lSO2emis + mfox*mSO2emis + hfox*hSO2emis)/ &
         (lSO2emis+mSO2emis+hSO2emis)

    if (lnuc.eq.0.d0 .and. mnuc.eq.0.d0 .and. hnuc.eq.0.d0) then
       ! no nucleation
       omasspp = 0.d0
       omdiam = 0.d0
       onnew = 0.d0
       ofnew = 0.d0
    else
       ! overall mass per particle is average weighted by number of 
       ! particles, disregarding no nucleation cases
       omasspp = (lmasspp*lSO2emis*lnnew*lnuc + &
                  mmasspp*mSO2emis*mnnew*mnuc + &
                  hmasspp*hSO2emis*hnnew*hnuc)/ &
                  (lSO2emis*lnnew*lnuc+mSO2emis*mnnew*mnuc+hSO2emis*hnnew*hnuc)
       
       ! overall nnew is total new particles over total emissions
       onnew = (lnnew*lSO2emis + mnnew*mSO2emis + hnnew*hSO2emis)/ &
            (lSO2emis+mSO2emis+hSO2emis)

       ! overall fnew is from closure
       ! note that fox*SO2emis weighted average would give the same answer
       tmpmasspp = omasspp
       tmpnnew = onnew
       call getfnew(ofox,tmpmasspp,tmpnnew,omasspp,onnew,ofnew)     

       ! overall mdiam just calculated from omasspp
       omdiam = getmdiam(omasspp)

    endif

  end subroutine P6_splitemis


  double precision function getfox(d,NOxemis,dswrf,vg,blh,bgNOx,ifoxp)
    ! returns fox, fraction of the emitted SO2 the has been oxidized at
    ! the defined distance from the source.

    !-----INPUTS----------------------------------------------------------

    double precision d       ! distance from source [m]
    double precision NOxemis ! NOx emissions [kg N/s]
    double precision dswrf   ! downward shortwave radiative flux [W/m]
    double precision vg      ! mean boundary layer wind speed [m/s]
    double precision blh     ! boundary layer height [m]
    double precision bgNOx   ! background NOx concentration [ppb]
    double precision, optional:: ifoxp(4) ! fitted parameters
    ! ifoxp should never be given as input, except by the other functions
    ! in this same module, where an effective fox value is calculated.

    !-----VARIABLE DECLARATIONS-------------------------------------------
    
    double precision time   ! time since emission [s]
    double precision effNOx ! effective in-plume NOx concentration [ppb]
    double precision effOH  ! effective OH concentration [molec cm-3]
    double precision p1 ! first polynomial; estimates shape of NOx vs. OH
    double precision p2 ! second polynomial; scales OH based on solar zenith angle
    double precision x
    double precision y

    double precision foxp(4) ! fox fitted parameters

    !-----ADJUSTABLE PARAMETERS-------------------------------------------
    
    double precision, parameter:: S0 = 1370.d0 ! Solar constant at top of atmosphere [W m-2]
    double precision, parameter:: Tr = 0.76d0  ! Transmittance through atmosphere

    double precision dfoxp(4) ! fox fitted parameters   
    double precision p1p(7)   ! parameters pertaining to the primary polynomial
    double precision p2p(4)   ! parameters for the second polynomial
   
    data dfoxp / -1.64966180d-10,   7.90402597d-01,   7.72321067d-01,   1.44390208d-08 /
    data p1p /-0.014d0, 0.0027d0, 0.1713d0, -0.0466d0, -0.7893d0, -0.1739d0, 6.9414d0/
    data p2p /-1345d0, 4002.d0, -471.8d0, 42.72d0/

    !-----CODE--------------------------------------------------------------

    if (present(ifoxp)) then
       foxp = ifoxp
    else
       foxp = dfoxp
    endif

    time = d/vg

    effNOx=bgNOx+foxp(4)*NOxemis*vg**noxfit(1)*blh**noxfit(2)*(time)**noxfit(3)

    ! The parameterization used here to calculate OHeff is described in
    ! appendix A of:
    ! Stevens, R. G., et al.,
    ! Nucleation and growth of sulfate aerosol in coal-fired power plant plumes:
    ! sensitivity to background aerosol and meteorology,
    ! Atmospheric Chemistry and Physics, doi:10.5194/acp-12-189-2012, 2012.

    x = log10(effNOx)-0.195
    y = dswrf/(S0*Tr)

    p1 = p1p(1)*x**6+p1p(2)*x**5+p1p(3)*x**4+p1p(4)*x**3+p1p(5)*x**2+p1p(6)*x+p1p(7)
    p2 = (p2p(1)*y**3+p2p(2)*y**2+p2p(3)*y+p2p(4))*1.E4
    
    effOH = (0.82 * 10.**(p1*log10(p2)/6.8))
   
    getfox = 1.-exp(foxp(1)*effOH**foxp(2)*time**foxp(3))

    return
  end function getfox


  logical function getnuc(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)
    ! returns True/False whether we expect to see not-insignificant
    ! new-particle formation.

    !-----INPUTS----------------------------------------------------------

    double precision d       ! distance from source [m]
    double precision SO2emis ! SO2 emissions [kg/s]
    double precision NOxemis ! NOx emissions [kg N/s]
    double precision CS      ! background condensation sink [s-1]
    double precision DSWRF   ! downward shortwave radiative flux [W/m]
    double precision vg      ! mean boundary layer wind speed [m/s]
    double precision BLH     ! boundary layer height [m]
    double precision bgSO2   ! background SO2 concentration [ppb]
    double precision bgNOx   ! background NOx concentration [ppb]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    double precision time     ! time since emission [s]
    double precision effSO2   ! effective in-plume SO2 concentration [ppb]
    double precision effNOx   ! effective in-plume NOx concentration [ppb]
    double precision nucp     ! nucleation predictor
    
    !-----ADJUSTABLE PARAMETERS-------------------------------------------
    
    double precision, parameter:: nucp_cutoff = 2.98841470581d+14
    double precision p(6)    ! list of fitted parameters used here

    data p / 4.35d0, 1.92d0, 3.28d0, -1.24d0, -3.48, 5.64d0/

    !-----CODE--------------------------------------------------------------

    time = d/vg

    effSO2 = bgSO2 + &
         10.d0**p(1)*SO2emis*vg**so2fit(1)*BLH**so2fit(2)*(time)**so2fit(3)

    effNOx = bgNOx + &
         10.d0**p(6)*NOxemis*vg**noxfit(1)*BLH**noxfit(2)*(time)**noxfit(3)

    ! avoid divide by zero:
    if (CS < 1e-5) then ! less than 1% of typical clean marine conditions
       getnuc = .True.
       ! unless someone chooses to call this directly, we do not have
       ! to worry about zero effNOx, because effNOx>=bgNOx>=minbgNOx
    else
       nucp = effSO2**p(2)*DSWRF**p(3)*effNOx**p(4)*CS**p(5)
       getnuc = (nucp .gt. nucp_cutoff)
    endif

    return
  end function getnuc


  double precision function getmasspp(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)

    ! returns masspp, the mean mass per particle for the newly formed particles

    !-----INPUTS----------------------------------------------------------

    double precision d       ! distance from source [m]
    double precision SO2emis ! SO2 emissions [kg/s]
    double precision NOxemis ! NOx emissions [kg N/s]
    double precision CS      ! background condensation sink [s-1]
    double precision DSWRF   ! downward shortwave radiative flux [W/m]
    double precision vg      ! mean boundary layer wind speed [m/s]
    double precision BLH     ! boundary layer height [m]
    double precision bgSO2   ! background SO2 concentration [ppb]
    double precision bgNOx   ! background NOx concentration [ppb]

    !-----VARIABLE DECLARATIONS-------------------------------------------
    
    double precision locfox   ! local effective fraction of SO2 oxidized [unitless]
    double precision time     ! time since emission [s]
    double precision effSO2   ! effective in-plume SO2 concentration [ppb]

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    double precision foxp(4) ! list of fitted parameters passed to fox
    double precision p(7)    ! list of fitted parameters used here

    data foxp /-1.29652905d-06, 6.92474330d-01, 2.92853444d-01, 2.13849343d+07/
    data p / 1.47496900d-27, 1.51723205d+00, 1.09357728d+00,-6.17290992d-01, &
             9.68490330d-01, 2.60502969d+06, 4.07112024d-23/

    !-----CODE--------------------------------------------------------------

    time = d/vg

    locfox = getfox(d,NOxemis,DSWRF,vg,BLH,bgNOx,foxp)
    
    effSO2 = bgSO2 + p(6)*SO2emis*vg**so2fit(1)*BLH**so2fit(2)*time**so2fit(3)

    getmasspp = p(1)*locfox**p(2)*effSO2**p(3)*CS**p(4)*time**p(5)+p(7)

    return
  end function getmasspp


  double precision function getmdiam(masspp)
    ! returns the particle number-median diameter in microns based on the
    ! mass per particle

    !-----INPUTS----------------------------------------------------------
    
    double precision masspp ! mass per particle of newly-formed
                            ! sulphate particles [kg]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    double precision mass_meand ! mass mean diameter

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    double precision, parameter:: dens = 1770.d0 
    ! density [kg/m3], assuming (NH4)2SO4
    double precision, parameter:: pi = 3.1415926535897931d0
    double precision, parameter:: sig = 1.4d0 
    ! geometric std. deviation [unitless]
    ! we assume a single lognormal mode, here

    !-----CODE--------------------------------------------------------------

    mass_meand = 1e6*(masspp/dens*6.d0/pi)**(1.d0/3.d0)

    getmdiam = mass_meand*exp(-1.5d0*log(sig)**2.d0)

    return
  end function getmdiam


  double precision function getnnew(d,SO2emis,NOxemis,CS,DSWRF,vg,BLH,bgSO2,bgNOx)
    ! returns nnew, the number of new particles per kg SO2 emitted

    !-----INPUTS----------------------------------------------------------

    double precision d       ! distance from source [m]
    double precision SO2emis ! SO2 emissions [kg/s]
    double precision NOxemis ! NOx emissions [kg N/s]
    double precision CS      ! background condensation sink [s-1]
    double precision DSWRF   ! downward shortwave radiative flux [W/m]
    double precision vg      ! mean boundary layer wind speed [m/s]
    double precision BLH     ! boundary layer height [m]
    double precision bgSO2   ! background SO2 concentration [ppb]
    double precision bgNOx   ! background NOx concentration [ppb]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    double precision locfox ! local effective fraction of SO2 oxidized [unitless]
    double precision time   ! time since emission [s]

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    double precision foxp(4) ! list of fitted parameters passed to fox
    double precision p(7)    ! list of fitted parameters used here

    data foxp /-3.54855422d-15, 7.13304235d-01, 1.93747558d+00, 1.24321647d+06/
    data p / 6.93853928d+23, 9.94909098d-01, 2.49960504d-01,-1.27968905d-01, &
            -4.41706268d+00, 1.44126017d-01, 1.73637370d-01/

    !-----CODE--------------------------------------------------------------

    time = d/vg
    
    locfox = getfox(d,NOxemis,DSWRF,vg,BLH,bgNOx,foxp)

    getnnew = p(1)*locfox**p(2)*bgSO2**p(3)*SO2emis**p(4)*exp(p(5)*CS**p(6)*(time)**p(7))+1.

    return
  end function getnnew


  subroutine getfnew(fox,imasspp,innew,omasspp,onnew,fnew)
    ! returns fnew, the fraction of formed H2SO4 mass in new particles 
    ! (instead of pre-existing particles), and adjusts masspp and nnew
    ! for closure if necessary.

    !-----INPUTS----------------------------------------------------------

    ! initial values of the following:
    double precision, intent(in):: fox
    ! fraction of emitted SO2 oxidized [unitless]
    double precision, intent(in):: imasspp
    ! mass per particle of newly-formed sulphate particles [kg]
    double precision, intent(in):: innew
    ! number of new particles per kg SO2 emitted [kg-1 SO2]

    !-----OUTPUTS---------------------------------------------------------

    ! output values of the following:
    double precision, intent(out):: omasspp
    ! mass per particle of newly-formed sulphate particles [kg]
    double precision, intent(out):: onnew
    ! number of new particles per kg SO2 emitted [kg-1 SO2]
    double precision, intent(out):: fnew
    ! fraction of H2SO4 formed that forms new particles [unitless]

    !-----VARIABLE DECLARATIONS-------------------------------------------

    double precision minmasspp ! minimum value allowed for masspp

    !-----ADJUSTABLE PARAMETERS-------------------------------------------

    double precision, parameter:: m_SO2 = 64.066d-3 ! molar mass of SO4 [kg/mol]
    double precision, parameter:: m_H2SO4 = 98.08d-3 ! molar mass of H2SO4 [kg/mol]
    double precision, parameter:: N_Av = 6.02214129d23 ! Avagadro's constant [1/mol]

    !-----CODE--------------------------------------------------------------

    minmasspp = 2.*m_H2SO4/N_Av ! mass of two H2SO4 molecules
    
    omasspp = imasspp
    onnew = innew

    fnew = (omasspp*onnew)/fox*m_SO2/m_H2SO4

    if (fnew > 1.) then
       ! where cfnew > 1., reduce it to one, and split the difference
       ! among masspp and nnew for closure.
       omasspp = omasspp/fnew**(0.5)
       onnew   = onnew/fnew**(0.5)

       ! I don't think this condition is possible, but just in case...
       if (omasspp < minmasspp) then
          onnew = onnew*omasspp/minmasspp
          omasspp = minmasspp
       endif

       fnew = 1.
    endif
     
    if (abs((omasspp*onnew)/(fox*fnew)*m_SO2/m_H2SO4 -1.)>eps) then
       print*, 'closure not achieved in get_fnew'
       print*, 'fnew',fnew
       print*, 'current ratio',(omasspp*onnew)/(fox*fnew)*m_SO2/m_H2SO4
    endif

    return
  end subroutine getfnew

end module P6mod
