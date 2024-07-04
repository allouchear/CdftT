include $(LIBCDFTTDIR)/CONFIG


BASIS = Basis
BECKE = Becke
CUBE = Cube
ORBITALS = Orbitals
UTILS = Utils
COMMON  = Common
JOBCONTROL  = JobControl

BASISOBJ  = $(BASIS)/*.o
BECKEOBJ  = $(BECKE)/*.o
ORBITALSOBJ  = $(ORBITALS)/*.o
UTILSOBJ  = $(UTILS)/*.o
CUBEOBJ  = $(CUBE)/*.o
COMMONSOBJ  = $(COMMON)/*.o
JOBCONTROLOBJ  = $(JOBCONTROL)/*.o

SUBDIRS = \
	$(COMMON) $(CUBE) $(BASIS) $(BECKE) $(ORBITALS) $(UTILS) $(JOBCONTROL)

OBJECTS = \
	$(COMMONOBJ) $(CUBEOBJ) $(BASISOBJ) $(BECKEOBJ) $(ORBITALSOBJ) $(UTILSOBJ) $(JOBCONTROLOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemilib

cchemilib: $(SUBDIRS)
	$(CC) -fPIC -shared $(LDFLAGS) $(OBJECTS) -o $(LIBCDFTTDIR)/lib/libcdftt.so  $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
	@$(RM) p
	@$(RM) pp
	@$(RM) p.cpp
	@$(RM) *.log
	@$(RM) *.out
	@$(RM) *.txt
