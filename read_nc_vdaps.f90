subroutine read_nc_vdaps(file_nc,ndiff,nvl,nt_nc,dt_nc,lats_nc,lons_nc,x,y,values_nc)
USE netcdf

type :: NC_DATA
 real*8, allocatable :: lat(:),lon(:),pre(:,:,:)
end type NC_DATA

type(NC_DATA) :: data_nc
real*8,dimension(NMAX) :: lats_nc,lons_nc
real*8,dimension(NMAX,nt_nc) :: values_nc
real*4,dimension(NMAX) :: x,y
real*4        :: lats_temp,lons_temp,x_temp,y_temp,data_temp
integer       :: xtype,ndim,nDims,nvar,varid
integer       :: dimids(3)
character*50  :: latname,lonname,timname,vname
character*150 :: file_nc

COMMON /NPAR/ NROW,NCOL,NMAX,CL,NCL
COMMON /LAMCON/ R,SLAT1,SLAT2,SN,SF,RO,OLAT,OLON,XO,YO,IFIRST,DD
COMMON /CALCUL/ PI,DEGRAD,RADDEG
COMMON /CHARDIR/ FCST_DIR,VDAPS_DIR,LDAPS_DIR,RDAPS_DIR,INTPOL,MGF

!stat=nf90_open(trim(file_nc),nf90_nowrite,ncid)
!write(*,*) stat

call check(nf90_open(trim(file_nc),nf90_nowrite,ncid))
call check(nf90_inquire(ncid,ndim,nvar))
call check(nf90_inquire_dimension(ncid,1,latname,nlat))
call check(nf90_inquire_dimension(ncid,2,lonname,nlon))
call check(nf90_inquire_dimension(ncid,3,timname,ntim))

write(*,*) ntim,nlat,nlon

spval=-99.9

allocate(data_nc%lat(nlat),data_nc%lon(nlon),data_nc%pre(nlon,nlat,ntim))

do i=1,nvar
 call check(nf90_inquire_variable(ncid,i,vname,xtype,nDims,dimids))
 call check(nf90_inq_varid(ncid,vname,varid))
 if(vname == "latitude") then
  call check(nf90_get_var(ncid,varid,data_nc%lat))
 elseif(vname == "longitude") then
  call check(nf90_get_var(ncid,varid,data_nc%lon))
 elseif(vname == "precipitation") then
  call check(nf90_get_var(ncid,varid,data_nc%pre))
 endif
 !write(*,*) trim(vname),varid,nDims,dimids(1),dimids(2),dimids(3)
enddo
do it=1,nt_nc
 nvl=0
 ny_st=YO-NCL;ny_re=YO+NROW+NCL
 nx_st=XO-NCL;nx_re=XO+NCOL+NCL
 do ilat=1,nlat
  do ilon=1,nlon
   lats_temp=data_nc%lat(ilat);lons_temp=data_nc%lon(ilon)
   data_temp=data_nc%pre(ilon,ilat,it)
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
   !endif
  enddo
 enddo
enddo

!write(*,*) nt_nc,nlat,nlon
!do it=1,1 !nt_nc
!do inv=1,nvl
! write(*,100) x(inv),y(inv),values_nc(inv,it)
!enddo
!enddo
!do inv=1,nvl
! temp_sum=0
! do it=49,54
!  temp=values_nc(inv,it)
!  if(temp.ge.0) then
!   temp_sum=temp_sum+temp
!  endif
! enddo
! write(*,100) x(inv),y(inv),temp_sum
!enddo

100 format(1000f10.2)
return
end
