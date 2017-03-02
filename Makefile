F90=ifort
FCFLAGS=-O3 -qopenmp

TARGET= SEM_CHANNEL
OBJECT= SEM_module.f90 SEM_main.f90 SEM_read.f90 SEM_write.f90

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean :
	rm -f *.o
	rm SEM_CHANNEL sem_module.mod
