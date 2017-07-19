COMPILER = gfortran
FLAGS = -O3 # this is an optimizer which makes loops more efficient
EXECUTABLE = 1d
# Object Files for build
OBJS = \
constants.o \
variables.o 
$(EXECUTABLE) : oned_euler.f90 $(OBJS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) oned_euler.f90 $(OBJS)
# Object dependencies and compilation
constants.o : constants.f90
	$(COMPILER) $(FLAGS) -c constants.f90 
variables.o : variables.f90 constants.o
	$(COMPILER) $(FLAGS) -c variables.f90 

all: run
run: 1d
	./1d
.PHONY: all run
