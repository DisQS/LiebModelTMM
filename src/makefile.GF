#-----------------------------------------------------------------
#
# tmseXdXX
#
#-----------------------------------------------------------------
#
# $Header: /home/cvs/phsht/tmseXd/src/makefile.GF,v 1.1 2016/10/02 19:13:05 phsht Exp $
#
# $Log: makefile.GF,v $
# Revision 1.1  2016/10/02 19:13:05  phsht
# now gfortran version
#
# Revision 1.1.1.1  2012/09/07 10:00:09  phsht
# tmseXd
#
# Revision 1.1  2011/05/06 08:13:09  phsht
# 1st installement
#
# Revision 1.8  2010/11/11 11:16:25  phrkaj
# Renamed files to get rid of any mention of SB.
#
# Revision 1.7  2010/11/04 15:14:02  phrkaj
# Added 3D stuff. Pre-debug
#
# Revision 1.6  2010/10/29 09:29:25  phrkaj
# *** empty log message ***
#
# Revision 1.5  2010/10/28 17:12:36  phrkaj
# *** empty log message ***
#
# Revision 1.4  2010/10/25 15:41:33  phsht
# small changes to remove a "malloc/glibc" error
#
# Revision 1.3  2010/10/25 11:11:28  phrkaj
# Back to pgi
#
# Revision 1.2  2010/10/24 14:18:50  phrkaj
# Changed compiler to gfortran for use at home. Added -fno-range-check to get it to work with gfortran.
#
# Revision 1.1.1.1  2010/10/22 12:23:38  phsht
# ParaTMM
#
#-----------------------------------------------------------------#

F77 = gfortran -ffixed-form -std=legacy
F90 = gfortran
#F90 = gfortran

#F77FLAGS = # For debugging
F77FLAGS = -fbounds-check -fno-range-check -Warray-bounds  # For optimizing
#F77FLAGS = # For profiling

F90FLAGS = -g -std=legacy -fbounds-check -fno-range-check -fimplicit-none -Warray-bounds -Wall # For debugging
#F90FLAGS = -ffree-line-length-none -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas # For optimizing
#F90FLAGS = -fno-range-check # for profling

CC = gcc
CCFLAGS = -O2

LIBDIR = $(HOME)/f77/lib
LIBFLAGS  =  -llapack -lblas

PRINT1 = a2ps -g -o print.ps
PRINT2 = > print.ps # | fleet -odouble

OBJ = common_modules.o application_modules.o random.o util.o main.o lieb21.o lieb22.o lieb23.o lieb31.o lieb32.o lieb33.o inout.o etime.o 

%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	tmseLMxD.GF
 
tmseLMxD.GF: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ) $(LIBFLAGS)

print:
	$(PRINT1) *.f90

clean:
	rm -f core *.mod *.o





