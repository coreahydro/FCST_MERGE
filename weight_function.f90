subroutine weight_function(nt,wNWP,wRDR)
real*4,dimension(0:nt) :: wNWP,wRDR
character*150 :: FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG
COMMON /CHARDIR/ FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

!TP=36.0 !! The time when (wNWP == 1) and (wRDR == 0)
!ALPHA_HTN=-0.1;BETA_HTN=1.1;GAMMA_HTN=1.0;DELTA_HTN=3.0
!ALPHA_SIG=16.0
!ALPHA_EQU=0.5
TC=TP/2 !! Time to cross between wNWP and wRDR
do it=0,nt
 IF(MGF.eq."SCN") then
!  t=2*PI/(nt*4)*it !! Period of function is nt
  if(it.le.TP) then
   t=2*PI/(TP*4)*it !! Period of function is tp
   wNWP(it)=(sin(t))**2
   wRDR(it)=(cos(t))**2
  else
   wNWP(it)=1
   wRDR(it)=0
  endif
 elseif(MGF.eq."HTN") then
  if(it.le.TP) then
   t=2*PI/(TP*8)*(it-TC)
   wNWP(it)=ALPHA_HTN+(BETA_HTN-ALPHA_HTN)/2.0*(1+tanh(GAMMA_HTN*DELTA_HTN*t))
   wRDR(it)=1-wNWP(it)
  else
   wNWP(it)=1
   wRDR(it)=0
  endif
 elseif(MGF.eq."SIG") then
  if(it.le.TP) then
   t=2*PI/(TP*8)*(it-TC)
   wNWP(it)=1/(1+exp(-1*ALPHA_SIG*t))
   wRDR(it)=1-wNWP(it)
  else
   wNWP(it)=1
   wRDR(it)=0
  endif
 elseif(MGF.eq."SMP") then
  if(it.le.TC) then
   wNWP(it)=0.0
   wRDR(it)=1.0
  else
   wNWP(it)=1.0
   wRDR(it)=0.0
  endif
 elseif(MGF.eq."EQU") then
  if(it.le.TP) then
   wNWP(it)=ALPHA_EQU
   wRDR(it)=1-ALPHA_EQU
  else
   wNWP(it)=1
   wRDR(it)=0
  endif   
 endif
!! write(*,*) it,wNWP(it),wRDR(it)
enddo

return
end
