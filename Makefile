F90=ifort
FCFLAGS=-O3 -qopenmp

TARGET= SEM_PROCESS
OBJECT= flow_module.o SEM_module.o main.o SEM_read.o SEM_write.o SEM_setup.o\
		SEM_eddy_setting.o SEM_fluctuation_gen.o SEM_combine_slice.o\
		SEM_convection.o SEM_statistics.o

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean :
	rm -f *.o
	rm SEM_PROCESS *.mod
