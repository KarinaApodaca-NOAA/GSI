.SUFFIXES: .f90 .o

OBJ =  da_recursive_filters.o da_reporting.o da_tools_serial.o da_etkf.o da_gen_be.o da_control.o \
module_ffts.o da_be_spectral.o da_fft_initialize.o da_change_wind_variables.o \
da_horizontal_covariances.o da_cloud_fraction.o configure.o variable_types.o io_input_model.o \
io_output.o io_input.o io_output_application.o read_benc.o

AR = ar
ARFLAGS = cr
TARGET = libtools.a

library: $(OBJ)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJ)
	$(RANLIB) $(TARGET)
	
da_recursive_filters.o : da_recursive_filters.f90
da_control.o : da_control.f90
module_ffts.o : module_ffts.f90
da_reporting.o : da_reporting.f90 da_control.o
da_tools_serial.o : da_tools_serial.f90 da_control.o da_reporting.o
da_etkf.o : da_etkf.f90 da_control.o
da_gen_be.o : da_gen_be.f90 da_control.o da_reporting.o da_tools_serial.o 
da_be_spectral.o : da_be_spectral.f90 da_control.o da_reporting.o da_tools_serial.o
da_fft_initialize.o : da_fft_initialize.f90 da_control.o da_tools_serial.o
da_change_wind_variables.o : da_change_wind_variables.f90 module_ffts.o
da_horizontal_covariances.o : da_horizontal_covariances.f90
da_cloud_fraction.o : da_cloud_fraction.f90
configure.o : configure.f90
variable_types.o : variable_types.f90
io_input_model.o : io_input_model.f90 configure.o variable_types.o da_gen_be.o 
io_output.o : io_output.f90 configure.o variable_types.o
io_output_application.o : io_output_application.f90 io_input_model.o io_input.o io_output.o
io_input.o : io_input.f90 configure.o variable_types.o
read_benc.o : read_benc.f90 

.f90.o :
	rm -f  $@
	$(CPP) $(CPPFLAGS) $(OMPCPP) $*.f90  > $*.f
	$(FC) -o $@ -c -I../blas -I../lapack -I$(NETCDFPATH)/include $(FCFLAGS) $*.f

clean:
	rm -f $(OBJ) $(TARGET) *.f *.mod

superclean: clean
