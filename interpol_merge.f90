subroutine interpol_merge(nt1,z1,nr2,nt2,x2,y2,z2,wNWP,wRDR,z_merge)
parameter(num_wt=130)
real*4,dimension(NMAX) :: x2,y2
real*8,dimension(NMAX,0:nt1) :: z1
real*8,dimension(NMAX,0:nt2) :: z_merge
real*8,dimension(NMAX,nt2) :: z2,z2_itp
real*4,dimension(NMAX,num_wt) :: wt
real*4,dimension(0:nt2) :: wNWP,wRDR
integer,dimension(NMAX,num_wt) :: n_z2
integer,dimension(NMAX) :: n_CL
integer,dimension(NMAX) :: near_num
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG,start_time
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

!!do ir=1,nr2
!! write(*,*) ir,x2(ir),y2(ir),z2(ir,1)
!!enddo

!!call cpu_time(end_time1)
!!cputime1=end_time1-start_time
!!write(*,*) "Time : ",cputime1/60," min"

 call read_merge_par(num_wt,near_num,n_z2,wt,n_CL)

!! IF want to calculate parameters for merging and write the parameters, close upper subroutine and open below subroutine
!! But the time to process will be so long (about 4 min) 
!! call cal_merge_par(num_wt,nr2,x2,y2,near_num,n_z2,wt,n_CL) 

!!write(*,*) INTPOL,MGF

!!call cpu_time(end_time2)
!!cputime2=end_time2-start_time
!!write(*,*) "Time : ",cputime2/60," min"

IF(INTPOL.eq."IDW") then
 do it=0,nt2
  !!write(*,*) it,wRDR(it),wNWP(it)
  if(it.eq.0) then ! it == 0, it(0) is observation
   do ixy=1,NMAX
    z1_temp=z1(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)
   enddo

  elseif(it.gt.0.and.it.le.nt1) then ! it < 0 && it <= nt1
   do ixy=1,NMAX
    nr_temp=n_CL(ixy)
    sum_wt=0;temp_z=0
    do ir=1,nr_temp
     temp_wt=wt(ixy,ir)
     ntemp=n_z2(ixy,ir)
     sum_wt=sum_wt+temp_wt
     temp_z=temp_z+temp_wt*z2(ntemp,it)
!     if(it.eq.1.and.ixy.eq.575031) write(*,*) ir,temp_wt,sum_wt,ntemp,z2(ntemp,it),temp_z
    enddo

    if(sum_wt.gt.0) then
     z2_itp(ixy,it)=temp_z/sum_wt
    else
     z2_itp(ixy,it)=0
    endif
      
    z1_temp=z1(ixy,it)
    z2_itp_temp=z2_itp(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)+z2_itp_temp*wNWP(it)
!    if(it.eq.12) write(*,100) it,wRDR(it),wNWP(it),z1_temp,z2_itp_temp,z_merge(ixy,it)
   enddo

  elseif(it.gt.nt1.and.it.le.nt2) then ! it < nt1 && it <= nt2
   do ixy=1,NMAX
    z2_itp_temp=z2_itp(ixy,it)
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z2_itp_temp*wNWP(it)
    !if(it.eq.37) write(*,100) it,wRDR(it),wNWP(it),z1_temp,z2_itp_temp,z_merge(ixy,it)
   enddo

  endif

 enddo
 
elseif(INTPOL.eq."SMT") then
 do it=0,nt2
  if(it.eq.0) then ! it == 0, it(0) is observation
   do ixy=1,NMAX
    z1_temp=z1(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)
   enddo 	
   
  elseif(it.gt.0.and.it.le.nt1) then ! it < 0 && it <= nt1
   do ixy=1,NMAX
    nr_temp=n_CL(ixy)
    sum_wt=0;temp_z=0
    do ir=1,nr_temp
     temp_wt=1
     ntemp=n_z2(ixy,ir)
     sum_wt=sum_wt+temp_wt
     temp_z=temp_z+temp_wt*z2(ntemp,it)
     !if(it.eq.1.and.ixy.eq.575031) write(*,*) ir,temp_wt,sum_wt,ntemp,z2(ntemp,it),temp_z
    enddo

    if(sum_wt.gt.0) then
     z2_itp(ixy,it)=temp_z/sum_wt
    else
     z2_itp(ixy,it)=0
    endif
    
    z1_temp=z1(ixy,it)
    z2_itp_temp=z2_itp(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)+z2_itp_temp*wNWP(it)
    !if(it.eq.18) write(*,100) it,wRDR(it),wNWP(it),z1_temp,z2_itp_temp,z_merge(ixy,it)      
   enddo
   
  elseif(it.gt.nt1.and.it.le.nt2) then ! it < nt1 && it <= nt2
   do ixy=1,NMAX
    z2_itp_temp=z2_itp(ixy,it)
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z2_itp_temp*wNWP(it)
   enddo
    
  endif
   
 enddo
 
else if(INTPOL.eq."NST") then
 do it=0,nt2
  !write(*,*) it
  if(it.eq.0) then ! it == 0, it(0) is observation
   do ixy=1,NMAX
    z1_temp=z1(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)
   enddo  
  
  elseif(it.gt.0.and.it.le.nt1) then ! it < 0 && it <= nt1  
   do ixy=1,NMAX
    nr_temp=near_num(ixy)     
    z2_itp(ixy,it)=z2(nr_temp,it)
    z1_temp=z1(ixy,it)
    z2_itp_temp=z2_itp(ixy,it)
    if(z1_temp.le.0) z1_temp=0
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z1_temp*wRDR(it)+z2_itp_temp*wNWP(it)
    !if(it.eq.18) write(*,100) it,wRDR(it),wNWP(it),z1(ixy,it),z2_itp(ixy,it),z_merge(ixy,it)
   enddo

  elseif(it.gt.nt1.and.it.le.nt2) then ! it < nt1 && it <= nt2
   do ixy=1,NMAX
    nr_temp=near_num(ixy) 
    z2_itp(ixy,it)=z2(nr_temp,it)
    z2_itp_temp=z2_itp(ixy,it)
    if(z2_itp_temp.le.0) z2_itp_temp=0
    z_merge(ixy,it)=z2_itp_temp*wNWP(it)
   enddo

  endif   
   
 enddo
endif
100 format(i5,200f10.2) 

!call cpu_time(end_time3)
!cputime3=end_time3-start_time
!write(*,*) "Time : ",cputime3/60," min":

return
end

subroutine read_merge_par(num_wt,near_num,n_z2,wt,n_CL)
real*4,dimension(NMAX,num_wt) :: wt
integer,dimension(NMAX,num_wt) :: n_z2
integer,dimension(NMAX) :: n_CL
integer,dimension(NMAX) :: near_num
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2

!!COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
!!COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
!!COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG,start_time
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

!write(CHAR_CL(1:2),'(i2.2)') int(CL)
open(1,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_NEAR_NUM_'//trim(CHAR_CL)//'KM.txt')
open(2,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_NUM_IN_'//trim(CHAR_CL)//'KM.txt')
open(3,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_WEIGT_IDW_'//trim(CHAR_CL)//'KM.txt')

nr_xy=0
do iy=1,NROW
 do ix=1,NCOL
  nr_xy=nr_xy+1
!  write(*,*) nr_xy
  read(1,*) ntempx,ntempy,near_num(nr_xy)
  read(2,*) ntempx,ntempy,n_temp,n_z2(nr_xy,1:n_temp)
  read(3,*) ntempx,ntempy,n_temp,wt(nr_xy,1:n_temp) 
  n_CL(nr_xy)=n_temp
 enddo
enddo 
close(1);close(2);close(3)

return
end

subroutine cal_merge_par(num_wt,nr2,x2,y2,near_num,n_z2,wt,n_CL)
real*4,dimension(NMAX) :: dist_near
real*4,dimension(NMAX) :: x2,y2
real*4,dimension(NMAX,num_wt) :: wt
integer,dimension(NMAX,num_wt) :: n_z2
integer,dimension(NMAX) :: n_CL
integer,dimension(NMAX) :: near_num
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG,start_time
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

!!write(CHAR_CL(1:2),'(i2.2)') int(CL)

open(1,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_NEAR_NUM_'//trim(CHAR_CL)//'KM.txt')
open(2,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_NUM_IN_'//trim(CHAR_CL)//'KM.txt')
open(3,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/'//trim(NWP_IN)//'_WEIGT_IDW_'//trim(CHAR_CL)//'KM.txt')

!!yc=YO
yc=0
nr_xy=0 ! Number of grid (Equal to NMAX)
do iy=1,NROW
 !!xc=XO
 xc=0
 do ix=1,NCOL
  n_temp=0 ! Number of grid in correlation length
  nr_xy=nr_xy+1
  dist_min=1000000.0
  do ir=1,nr2
   dist=sqrt(((xc-x2(ir))*DD)**2+((yc-y2(ir))*DD)**2)
   if(dist.le.CL) then
    if(dist.le.dist_min) then
     dist_min=dist
     near_num(nr_xy)=ir
     dist_near(nr_xy)=dist_min
    endif
    n_temp=n_temp+1
    temp_wt=1/(dist**2)
    if(dist.le.0.01) temp_wt=1000
    wt(nr_xy,n_temp)=temp_wt
    n_z2(nr_xy,n_temp)=ir
    !!if(ix.eq.620.and.iy.eq.210) write(*,*) ix,iy,ir,xc,yc,x2(ir),y2(ir),dist,temp_wt,z2(ir,1)
   endif
  enddo
  n_CL(nr_xy)=n_temp
  write(1,997) ix,iy,near_num(nr_xy)
  write(2,998) ix,iy,n_CL(nr_xy),n_z2(nr_xy,1:n_temp)
  write(3,999) ix,iy,n_CL(nr_xy),wt(nr_xy,1:n_temp)
  xc=xc+1
 enddo
 yc=yc+1
enddo
997 format(3i8)
998 format(1000i8)
999 format(3i8,1000f8.2)
close(1);close(2);close(3)

return
end

