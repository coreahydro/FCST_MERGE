subroutine write_nc(file_nc,nmax,nx,ny,nt,values_nc)
USE netcdf

real*8,dimension(nmax,0:nt) :: values_nc
real*8 temp
integer*2 data_nc(ny,nx,nt+1)
integer*4 dimids(3) !! Must integer*4
character*150  file_nc

do it=0,nt
 nr=0
 do iy=1,ny
  do ix=1,nx
   nr=nr+1
   temp=values_nc(nr,it)
   if(temp.lt.0.or.temp.gt.10000) then
    data_nc(iy,ix,it+1)=-999
   else    	
    data_nc(iy,ix,it+1)=temp*10
   endif
!   if(it.eq.18) write(*,*) nr,iy,ix,values_nc(nr,it),data_nc(iy,ix,it+1)
  enddo
 enddo
enddo

call check(nf90_create(trim(file_nc),nf90_clobber,ncid))
call check(nf90_def_dim(ncid,"NX",nx,nx_dimid))
call check(nf90_def_dim(ncid,"NY",ny,ny_dimid))
call check(nf90_def_dim(ncid,"NT",nt+1,nt_dimid))
dimids=(/ny_dimid,nx_dimid,nt_dimid/)
call check(nf90_def_var(ncid,"FCST_VDAPS",NF90_INT2,dimids,ID_FCVD)) !!integer*4=NF90_INT4, integer*2=NF90_INT2
call check(nf90_enddef(ncid))
call check(nf90_put_var(ncid,ID_FCVD,data_nc))
call check(nf90_close(ncid))

return
end
