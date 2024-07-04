include $(LIBCDFTTDIR)/CONFIG

BASIS = Basis
BECKE = Becke
ORBITALS = Orbitals
UTILS = Utils

BASISOBJ = $(BASIS)/*.o
BECKEOBJ  = $(BECKE)/*.o
ORBITALSOBJ  = $(ORBITALS)/*.o
UTILSOBJ = $(UTILS)/*.o

SUBDIRS = \
	$(BASIS) $(BECKE) $(ORBITALS) $(UTILS)

OBJECTS = \
	$(BASISOBJ) $(BECKEOBJ) $(ORBITALSOBJ) $(UTILSOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemilib

cchemilib: $(SUBDIRS)
	$(CC) -fPIC -shared $(LDFLAGS) $(OBJECTS) -o ../../../lib/libcdftt.so  $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
	@$(RM) p
	@$(RM) pp
	@$(RM) p.cpp
	@$(RM) *.log
	@$(RM) *.out
	@$(RM) *.txt
