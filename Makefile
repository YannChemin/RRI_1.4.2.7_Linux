TARGET = 0_rri_1_4_2_7_Linux

OBJECT = RRI_Mod.o RRI_Mod2.o RRI_Mod_Dam.o RRI_Mod_Tecout.o RRI.o RRI_Bound.o RRI_Dam.o RRI_Div.o RRI_DT_Check.o RRI_Evp.o RRI_GW.o RRI_Infilt.o RRI_Read.o RRI_Riv.o RRI_RivSlo.o RRI_Section.o RRI_Slope.o RRI_Sub.o RRI_Tecout.o RRI_TSAS.o

#F90 = ifort
#F90 = nfort-3.4.2
F90 = gfortran

#FFLAGS = -O -qopenmp
#FFLAGS = -O3 -fopenmp
# OFFLOADING
FFLAGS = -O3 -fopenmp -foffload=-lgfortran

INCLDIR = .
#INCLDIR = 
.SUFFIXES : .o .f90
.f90.o:
	${F90} ${FFLAGS} -c $< -I${INCLDIR}
${TARGET}: ${OBJECT}
	${F90} -o $@ ${OBJECT} ${FFLAGS}

mod:
	$(F90) $(FFLAGS) -c RRI_Mod.f90
	$(F90) $(FFLAGS) -c RRI_Mod2.f90
	$(F90) $(FFLAGS) -c RRI_Mod_Dam.f90
	$(F90) $(FFLAGS) -c RRI_Mod_Tecout.f90

clean:

	rm -f ${TARGET} *.L *.mod *.o
