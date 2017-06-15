F90 = ftn
OPTS = -O3
LIB = $(DFFTPACK)
OBJS = li_com.f90 fft_wrapper.f90 fcnt.f90 li.f90

li:  $(OBJS)
	$(F90) $(OPTS) $(OBJS) -o li $(LIB)

clean:
	rm -f li *.o *.mod
