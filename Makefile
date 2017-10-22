F90=ifort
FCFLAGS=-O3 -qopenmp

TARGET= SEM_PROCESS
OBJECT= math_module.o IO_module.o flow_module.o SEM_module.o main.o SEM_write.o

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean :
	rm -f *.o
	rm SEM_PROCESS *.mod
