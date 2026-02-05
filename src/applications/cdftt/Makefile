#LIBCDFTTDIR =/usr/home/lgardre-local/Work/CdftT

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

OBJECTS = cdftt.o \
	$(COMMONOBJ) $(CUBEOBJ) $(BASISOBJ) $(BECKEOBJ) $(ORBITALSOBJ) $(UTILSOBJ) $(JOBCONTROLOBJ)

all: compiling  cdftt

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.inp *~ cdftt *.mk

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

cdftt: 
	$(CC) $(LDFLAGS) $(OBJECTS) -o cdftt $(OMPLIB) -g -lm $(ACCLIB)  


-include Dep.mk

