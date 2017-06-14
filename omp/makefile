F90 = gfortran
OPTS = -fopenmp  -Ofast
#OPTS = -Wall -fbacktrace -fcheck=all
#F90 = ifort
#OPTS = -qopenmp -O3
LIB = -L/usr/local/lib/dfftpack -ldfftpack
OBJS = li_com.f90 fft_wrapper.f90 fcnt.f90 li.f90

li:  $(OBJS)
	$(F90) $(OPTS) $(OBJS) -o li $(LIB)

clean:
	rm -f li *.o *.mod
