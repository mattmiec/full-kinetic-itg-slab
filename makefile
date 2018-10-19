#ubuntu
FC_ubuntu= mpifort
FLAGS_ubuntu = -Ofast
#FLAGS_ubuntu = -fcheck=all -fbacktrace
#INC_ubuntu = -I/usr/local/include
#LIB_ubuntu = -L/usr/local/lib -lfftw3
LIB_ubuntu = -L/usr/local/lib/dfftpack -ldfftpack

#edison/cori
FC_nersc = ftn
FLAGS_nersc = -fast
#FLAGS_nersc = -check all -traceback
LIB_nersc = $(DFFTPACK)

#rmacc summit
#FC_summit = mpifort
#FLAGS_summit_cpu = -fast
#FLAGS_summit_gpu = -fast -acc -ta=tesla:cc35 -Minfo=accel
#INC_summit = -I$(CURC_FFTW_INC)
#LIB_summit = -L$(CURC_FFTW_LIB)  -lfftw3

SRCS = li_com.f90 fft_wrapper.f90 fcnt.f90 li.f90

li_ubuntu: $(SRCS)
	$(FC_ubuntu) $(FLAGS_ubuntu) $(INC_ubuntu) $(SRCS) $(LIB_ubuntu) -o li_ubuntu

li_cori:  $(SRCS)
	$(FC_nersc) $(FLAGS_nersc) $(SRCS) $(LIB_nersc) -o li_cori

li_edison:  $(SRCS)
	$(FC_nersc) $(FLAGS_nersc) $(SRCS) $(LIB_nersc) -o li_edison

#li_summit_cpu:  $(SRCS)
#	$(FC_summit) $(FLAGS_summit_cpu) $(INC_summit) $(SRCS) $(LIB_summit) -o li_summit_cpu

#li_summit_gpu:  $(SRCS)
#	$(FC_summit) $(FLAGS_summit_gpu) $(INC_summit) $(SRCS) $(LIB_summit) -o li_summit_gpu

clean:
	rm -f li_ubuntu li_cori li_edison li_summit_cpu li_summit_gpu *.o *.mod
