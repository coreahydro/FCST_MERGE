subroutine read_par(file_par)
character*150 :: file_par
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL,NT_FCST,NT_NWP
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

open(1,file=trim(file_par))
!! Information for input data !!
read(1,*);read(1,*);read(1,*)
read(1,*) NROW;read(1,*) NCOL
read(1,*) NT_FCST
read(1,*) NT_NWP
read(1,100) RES_DIR
read(1,100) FCST_DIR
read(1,100) VDAPS_DIR
read(1,100) LDAPS_DIR
read(1,100) RDAPS_DIR
NT_FCST=NT_FCST-1 !! 0 is observation 
NMAX=NROW*NCOL

!! Information for input data !!
read(1,*);read(1,*);read(1,*)
read(1,101) OUTFORM
read(1,100) OUT_DIR
!! Map Information !!
read(1,*);read(1,*);read(1,*)
read(1,*) CL
read(1,*) R;read(1,*) SLAT1;read(1,*) SLAT2
read(1,*) OLAT;read(1,*) OLON
read(1,*) XO;read(1,*) YO
read(1,*) DD

write(CHAR_CL(1:2),'(i2.2)') int(CL)
NCL=CL/DD
PI = ASIN(1.)*2.
DEGRAD = PI/180.
RADDEG = 180./PI
SLAT1 = SLAT1*DEGRAD
SLAT2 = SLAT2*DEGRAD
OLAT = OLAT*DEGRAD
OLON = OLON*DEGRAD
R = R/DD

!! Information for Merge Function !!
read(1,*);read(1,*);read(1,*)
read(1,*) TP
read(1,*) ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN
read(1,*) ALPHA_SIG
read(1,*) ALPHA_EQU

!write(*,*) trim(RES_DIR)
!write(*,*) trim(FCST_DIR)
!write(*,*) trim(VDAPS_DIR)
!write(*,*) trim(LDAPS_DIR)
!write(*,*) trim(RDAPS_DIR)
!write(*,*) trim(OUTFORM)
!write(*,*) trim(OUT_DIR)
!write(*,*) ALPHA_SIG
!write(*,*) ALPHA_EQU

100 format(a150)
101 format(a2)

return
end
