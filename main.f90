!! This program merges very short term precipitation forecasting data(FCST) with numerical weahter prediction data
!! Programed by Jungsoo YOON beloing to KICT
!! If you have any question, please send email to "jungsooyoon@kict.re.kr"
!!USE grib_api
USE netcdf
!! ny=419,nx=491 RDAPS
logical :: lexist1,lexist2
real*8,dimension(:,:),allocatable :: data_fcst,data_vdaps,data_klaps,data_merge
real*8,dimension(:),allocatable :: lats_vdaps,lons_vdaps
real*4,dimension(:),allocatable :: x_vdaps,y_vdaps
real*8,dimension(:),allocatable :: lats_klaps,lons_klaps
real*4,dimension(:),allocatable :: x_klaps,y_klaps
real*8,dimension(:),allocatable :: lats_ldaps,lons_ldaps,data_ldaps
real*4,dimension(:),allocatable :: x_ldaps,y_ldaps
real*8,dimension(:),allocatable :: lats_rdaps,lons_rdaps,data_rdaps
real*4,dimension(:),allocatable :: x_rdaps,y_rdaps
real*4,dimension(:),allocatable :: wNWP,wRDR
integer :: lagtime,diff_time
character*150                   :: file_par,file_fcst,file_vdaps,file_ldaps,file_rdaps,file_nwp
character*150                   :: file_merge,makedir
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2
character :: seltime*12,murge_func*3
character :: fcst_time*12,nwp_time*12,diff_time_char*3
character :: fcst_info*18,nwp_info*18

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL,NT_FCST,NT_NWP
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG,start_time
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

call cpu_time(start_time)
nget=1
call getarg(nget,file_par)
call read_par(file_par)

nget=nget+1
call getarg(nget,NWP_IN)

nget=nget+1
call getarg(nget,INTPOL)

nget=nget+1
call getarg(nget,MGF)

nget=nget+1
call getarg(nget,fcst_time)

nget=nget+1
call getarg(nget,nwp_time)

nget=nget+1
call getarg(nget,diff_time_char)
read(diff_time_char(1:3), '(i3.3)') diff_time

write(*,*) fcst_time," ",nwp_time,diff_time

nget=nget+1
call getarg(nget,file_fcst)

nget=nget+1
call getarg(nget,file_nwp)

!!write(*,*) INTPOL," ",MGF
!!write(*,*) trim(file_fcst)
!!write(*,*) trim(file_nwp)

allocate(wNWP(0:NT_NWP),wRDR(0:NT_NWP)) !0 is observation, so including 0
call weight_function(NT_NWP,wNWP,wRDR)

nget=nget+1
!!file_fcst=trim(FCST_DIR)//'/'//yyyy_k//mt_k//'/'//dd_k//'/RDR_QPF_FCST_'//kst//'.bin.gz'
!!file_fcst=trim(FCST_DIR)//'/RDR_QPF_FCST_'//kst//'.bin.gz'
!!write(*,*) trim(file_fcst)
inquire(file=trim(file_fcst),exist=lexist1) 
if(.not.lexist1) then
 write(*,*) " The fcst file does not exist!"
 write(*,*) " Stop the program"
 stop
else 
 write(*,*) "INPUT FILE1: ",file_fcst
endif
allocate(data_fcst(NMAX,0:NT_FCST)) !! 0 is observation 
call read_bin(file_fcst,NMAX,NT_FCST,data_fcst)
!!do it=0,0
!! do ir=1,NMAX
!!  write(*,*) it,ir,data_fcst(ir,it)
!! enddo
!!enddo

!!call LAMCPROJ(33.7,124.45,x_temp,y_temp,0,IERR)
!!xc_temp=x_temp+1050/2
!!yc_temp=y_temp+1050/2
!!call LAMCPROJ(temp_lat,temp_lon,xc_temp,yc_temp,1,IERR)
!!write(*,*) x_temp,y_temp
!!write(*,*) temp_lat,temp_lon

!! Make the directory of out data
!!makedir="mkdir -p "//trim(OUT_DIR)//'/'//yyyy_k//mt_k//'/'//dd_k
makedir="mkdir -p "//trim(OUT_DIR)
call system(makedir)

DT_VDAPS=10.0 ! Interval Time of VDAPS(minute)
NT_VDAPS=NT_NWP
DT_KLAPS=10.0 ! Interval Time of KLAPS(minute)
NT_KLAPS=NT_NWP
ndiff=diff_time/DT_VDAPS
if(ndiff.lt.0) then
 write(*,*) " VDPAS file is later than FCST file!"
 write(*,*) " Stop the program"
 stop
endif 
fcst_info="FCST "//fcst_time//" "
nwp_info=NWP_IN//" "//nwp_time

if(NWP_IN.eq.'VDAPS') then
 nget=nget+1
!! file_nwp=trim(VDAPS_DIR)//'/'//yyyy_u//mt_u//'/'//dd_u//'/vdps_10min_rain_'//utc//'.nc'
!! file_nwp=trim(VDAPS_DIR)//'/vdps_10min_rain_'//utc//'.nc'
! write(*,*) trim(file_vdaps)
 inquire(file=trim(file_nwp),exist=lexist2) 
 if(.not.lexist2) then
  write(*,*) " The VDPAS file does not exist!"
  write(*,*) " Stop the program"
  stop
 else 
  write(*,*) "INPUT FILE2: ",file_nwp
 endif
 allocate(lats_vdaps(NMAX),lons_vdaps(NMAX),x_vdaps(NMAX),y_vdaps(NMAX),data_vdaps(NMAX,NT_VDAPS))
 allocate(data_merge(NMAX,0:NT_VDAPS)) !! 0 is observation 
 call read_nc_vdaps(file_nwp,ndiff,nvl_vdaps,NT_VDAPS,DT_VDAPS,lats_vdaps,lons_vdaps,x_vdaps,y_vdaps,data_vdaps)
 call interpol_merge(NT_FCST,data_fcst,nvl_vdaps,NT_VDAPS,x_vdaps,y_vdaps,data_vdaps,wNWP,wRDR,data_merge) 
 if(OUTFORM.eq."BI") then
  !!file_merge=trim(OUT_DIR)//'/'//yyyy_k//mt_k//'/'//dd_k//'/MGRF_FCVD_'//INTPOL//'_'//MGF//'_'//fcst_time//'.bin.gz'
  file_merge=trim(OUT_DIR)//'/MGRF_FCVD_'//INTPOL//'_'//MGF//'_'//fcst_time//'.bin.gz'
  write(*,*) "OUTPUT FILE: ",file_merge
  call write_bin(file_merge,fcst_info,nwp_info,NMAX,NCOL,NROW,NT_VDAPS,data_merge,NT_FCST)
 elseif(OUTFORM.eq."NC") then
  !!file_merge=trim(OUT_DIR)//'/'//yyyy_k//mt_k//'/'//dd_k//'/MGRF_FCVD_'//INTPOL//'_'//MGF//'_'//fcst_time//'.nc'
  file_merge=trim(OUT_DIR)//'/MGRF_FCVD_'//INTPOL//'_'//MGF//'_'//fcst_time//'.nc'
  call write_nc(file_merge,NMAX,NCOL,NROW,NT_VDAPS,data_merge)
  call system("gzip "//trim(file_merge)) !! gzip netcdf out data
  write(*,*) "OUTPUT FILE: ",file_merge
 endif
elseif(NWP_IN.eq.'KLAPS') then
 nget=nget+1
 inquire(file=trim(file_nwp),exist=lexist2) 
 if(.not.lexist2) then
  write(*,*) " The KLAPS file does not exist!"
  write(*,*) " Stop the program"
  stop
 else 
  write(*,*) "INPUT FILE2: ",file_nwp
 endif
 allocate(lats_klaps(NMAX),lons_klaps(NMAX),x_klaps(NMAX),y_klaps(NMAX),data_klaps(NMAX,NT_KLAPS))
 allocate(data_merge(NMAX,0:NT_KLAPS)) !! 0 is observation 
 call read_nc_klaps(file_nwp,ndiff,nvl_klaps,NT_KLAPS,DT_KLAPS,lats_klaps,lons_klaps,x_klaps,y_klaps,data_klaps)
 call interpol_merge(NT_FCST,data_fcst,nvl_klaps,NT_KLAPS,x_klaps,y_klaps,data_klaps,wNWP,wRDR,data_merge) 
 if(OUTFORM.eq."BI") then
  file_merge=trim(OUT_DIR)//'/MGRF_FCKL_'//INTPOL//'_'//MGF//'_'//fcst_time//'.bin.gz'
  write(*,*) "OUTPUT FILE: ",file_merge
  call write_bin(file_merge,fcst_info,nwp_info,NMAX,NCOL,NROW,NT_KLAPS,data_merge,NT_FCST)
 elseif(OUTFORM.eq."NC") then
  file_merge=trim(OUT_DIR)//'/MGRF_FCKL_'//INTPOL//'_'//MGF//'_'//fcst_time//'.nc'
  call write_nc(file_merge,NMAX,NCOL,NROW,NT_KLAPS,data_merge)
  call system("gzip "//trim(file_merge)) !! gzip netcdf out data
  write(*,*) "OUTPUT FILE: ",file_merge
 endif

else
 write(*,*) "INPUT FILE2 is WRONG"
endif

call cpu_time(end_time)
cputime=end_time-start_time
!write(*,*) "Time : ",cputime/60," min"

stop
end

