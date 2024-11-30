subroutine read_nc_klaps(file_nc,ndiff,nvl,nt_nc,dt_nc,lats_nc,lons_nc,x,y,values_nc)
USE netcdf

type :: NC_DATA
 real*8, allocatable :: lat(:,:),lon(:,:),rainc(:,:,:),rainnc(:,:,:)
end type NC_DATA

type(NC_DATA) :: data_nc
real*8,dimension(NMAX) :: lats_nc,lons_nc
real*8,dimension(NMAX,nt_nc) :: values_nc
real*4,dimension(NMAX) :: x,y
real*4        :: lats_temp,lons_temp,x_temp,y_temp,data_temp
integer       :: xtype,ndim,nDims,nvar,varid,ntemp
integer       :: dimids(3)
character*50  :: latname,lonname,timname,datelen,vname
character*50  :: temp
character*150 :: file_nc
character*150 :: RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
character :: INTPOL*3,MGF*3,NWP_IN*5,CHAR_CL*2,OUTFORM*2

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL,NT_FCST,NT_NWP
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG,start_time
COMMON /CHARDIR/ RES_DIR,FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,OUT_DIR
COMMON /CHARPAR/ INTPOL,MGF,NWP_IN,CHAR_CL,OUTFORM
COMMON /FUNCPAR/ TP,ALPHA_HTN,BETA_HTN,GAMMA_HTN,DELTA_HTN,ALPHA_SIG,ALPHA_EQU

!stat=nf90_open(trim(file_nc),nf90_nowrite,ncid)
!write(*,*) stat

call check(nf90_open(trim(file_nc),nf90_nowrite,ncid))
call check(nf90_inquire(ncid,ndim,nvar))
call check(nf90_inquire_dimension(ncid,1,timname,ntim))
call check(nf90_inquire_dimension(ncid,2,datelen,ndln))
call check(nf90_inquire_dimension(ncid,3,lonname,nlon))
call check(nf90_inquire_dimension(ncid,4,latname,nlat))

!!call check(nf90_get_att(ncid,nf90_global,"WEST-EAST_GRID_DIMENSION", ntemp) )
!!call check(nf90_inq_attname(ncid,nf90_global,"title",ntemp))

spval=-99.9

allocate(data_nc%lat(nlat,nlon),data_nc%lon(nlat,nlon),data_nc%rainc(nlon,nlat,ntim),data_nc%rainnc(nlon,nlat,ntim))

open(1,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/latitude_klaps_kma.txt')
open(2,file=trim(RES_DIR)//'/'//trim(NWP_IN)//'/longitude_klaps_kma.txt')

read(1,*) nlat_max
read(2,*) nlon_max

if(nlat_max.ne.nlat) then
	write(*,*) " The number of south_north(latitude) is not equal to input information of latitude!"
	write(*,*) " Check file of 'latitude_klaps_kma.txt'"	
	write(*,*) " Stop the program"
	stop
endif

if(nlon_max.ne.nlon) then
	write(*,*) " The number of west_east(longitude) is not equal to input information of longitude!"
	write(*,*) " Check file of 'longitude_klaps_kma.txt'"	
	write(*,*) " Stop the program"
	stop
endif

do ilat=1,nlat
	read(1,*) data_nc%lat(ilat,1:nlon)
	read(2,*) data_nc%lon(ilat,1:nlon) !! Correction(2023.05.23.) 
enddo

close(1);close(2)

do i=1,nvar
 call check(nf90_inquire_variable(ncid,i,vname,xtype,nDims,dimids))
 call check(nf90_inq_varid(ncid,vname,varid))
 if(vname == "RAINC") then
  call check(nf90_get_var(ncid,varid,data_nc%rainc))
 elseif(vname == "RAINNC") then
  call check(nf90_get_var(ncid,varid,data_nc%rainnc))
 endif
! write(*,*) trim(vname),varid,nDims,dimids(1),dimids(2),dimids(3)
enddo

do it=1,nt_nc !! nt_nc(72) equal to ntim(73)-1 --> First rain data remove
 nvl=0
 ny_st=YO-NCL;ny_re=YO+NROW+NCL
 nx_st=XO-NCL;nx_re=XO+NCOL+NCL
!! write(*,*) nx_st,nx_re,ny_st,ny_re
 do ilat=1,nlat
  do ilon=1,nlon
   lats_temp=data_nc%lat(ilat,ilon);lons_temp=data_nc%lon(ilat,ilon)
   temp_rainc=data_nc%rainc(ilon,ilat,it+1)-data_nc%rainc(ilon,ilat,it)
   temp_rainnc=data_nc%rainnc(ilon,ilat,it+1)-data_nc%rainnc(ilon,ilat,it)
   data_temp=temp_rainc+temp_rainnc
   call LAMCPROJ(lats_temp,lons_temp,x_temp,y_temp,0,IERR)
   !!if(it.eq.1) write(*,*) lats_temp,lons_temp,x_temp,y_temp
   !!if (x_temp.ge.nx_st.and.x_temp.le.nx_re.and.y_temp.ge.ny_st.and.y_temp.le.ny_re) then
    nvl=nvl+1
    values_nc(nvl,it)=spval
    x(nvl)=x_temp;y(nvl)=y_temp
    lats_nc(nvl)=lats_temp;lons_nc(nvl)=lons_temp
!!    if(ntemp.ge.1) values_nc(nvl,it)=data_temp !*60.0/dt_nc ! Convert mm to mm/hr
    ntemp=it-ndiff
    if(ntemp.ge.1) values_nc(nvl,ntemp)=data_temp*60.0/dt_nc ! Convert mm to mm/hr
    !!if(it.eq.24) write(*,100) nvl,lats_nc(nvl),lons_nc(nvl),x(nvl),y(nvl),data_temp 
   !!endif
  enddo
 enddo
enddo

100 format(i5,1000f10.2)
return
end
