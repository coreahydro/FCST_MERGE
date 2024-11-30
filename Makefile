FF=gfortran
FC=gcc

NETCDFF_ROOT = /usr/local/netcdf-fortran
NETCDFC_ROOT = /usr/local/netcdf
#NETCDFF_ROOT = /usr/local/netcdf-fortran/4.4.3
#NETCDFC_ROOT = /usr/local/netcdf-c/4.4.0

INCLUDES    = -I${NETCDFF_ROOT}/include \
	      -I${NETCDFC_ROOT}/include \

LIBRARIES   = -L${NETCDFF_ROOT}/lib -lnetcdff -lz -lm \
	      -L${NETCDFC_ROOT}/lib -lnetcdf \

SRCS_C = read_bin_data.c\
         write_bin_data.c\

SRCS_F = main.f90\
     read_par.f90\
     read_nc_vdaps.f90\
     read_nc_klaps.f90\
     check.f90 \
     lamcproj.f\
     weight_function.f90\
     interpol_merge.f90\
     write_nc_data.f90\

DEPLOY_DIR      = "../BIN"

OBJS_C = ${SRCS_C:.c=.o}
OBJS_F = ${SRCS_F:.f=.o}

all : fcst_merge

.SUFFIXES : .o .c .f

fcst_merge : ${OBJS_C} ${OBJS_F}
	${FF} -fbounds-check -o $@ $^  $(INCLUDES) ${LIBRARIES}
	@cp -a $@ ${DEPLOY_DIR}

${OBJS_C}: ${SRCS_C} 
	${FC} -c  ${SRCS_C} $(INCLUDES) ${LIBRARIES}	
	
${OBJS_F}: ${SRCS_F}
	${FF} -c ${SRCS_F} $(INCLUDES) ${LIBRARIES}

clean :
	\rm -f *.c~ *.h~ *.o *.f~ Makefile~ fcst_merge
