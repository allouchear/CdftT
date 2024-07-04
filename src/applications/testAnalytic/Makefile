include $(LIBCDFTTDIR)/CONFIG

CFLAGS =  $(COMMONCFLAGS) 

BASIS = $(LIBCDFTTDIR)/src/lib/Basis
BECKE = $(LIBCDFTTDIR)/src/lib/Becke
CUBE = $(LIBCDFTTDIR)/src/lib/Cube
ORBITALS = $(LIBCDFTTDIR)/src/lib/Orbitals
UTILS = $(LIBCDFTTDIR)/src/lib/Utils
COMMON  = $(LIBCDFTTDIR)/src/lib/Common
JOBCONTROL  = $(LIBCDFTTDIR)/src/lib/JobControl

BASISOBJ  = $(BASIS)/*.o
BECKEOBJ  = $(BECKE)/*.o
ORBITALSOBJ  = $(ORBITALS)/*.o
UTILSOBJ  = $(UTILS)/*.o
CUBEOBJ  = $(CUBE)/*.o
COMMONOBJ  = $(COMMON)/*.o
JOBCONTROLOBJ  = $(JOBCONTROL)/*.o

OBJECTS = main.o \
	$(COMMONOBJ) $(CUBEOBJ) $(BASISOBJ) $(BECKEOBJ) $(ORBITALSOBJ) $(UTILSOBJ) $(JOBCONTROLOBJ)

all: compiling cdftt

%.o : %.cpp
	$(CC) -g -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) cdftt *o *.out *.inp *~ 

cdftt: 
	$(CC) $(LDFLAGS) $(OBJECTS) -o cdftt $(OMPLIB) -g -lm -lgomp 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk



-include Dep.mk

