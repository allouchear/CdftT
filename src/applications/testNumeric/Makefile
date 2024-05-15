#pour compiler, taper Taper ./compcdftt.sh

include $(LIBCDFTTDIR)/CONFIG

MAINDIR   = ./
ANALYTIC  = $(LIBCDFTTDIR)/src/analytic
NUMERIC   = $(LIBCDFTTDIR)/src/numeric
COMMON    = $(LIBCDFTTDIR)/src/common


ANALYTICOBJ = $(ANALYTIC)/*.o
NUMERICOBJ  = $(NUMERIC)/*.o
COMMONSOBJ  = $(COMMON)/*.o



MAINOBJ   = $(MAINDIR)/*.o
 
SUBDIRS = \
	$(MAINDIR) $(ANALYTIC) $(NUMERIC) $(COMMON) 

OBJECTS = \
	$(MAINOBJ) $(ANALYTICOBJ) $(ANALYTICOBJ) $(COMMONOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemi

cchemi: $(SUBDIRS)
	g++ $(LDFLAGS) $(OBJECTS) -o cchemi $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm -L $(LIBNNPDIR) -lnnp  -lgomp $(PYTHONLIB)
#	$(CC) $(LDFLAGS) $(OBJECTS) -o cchemi $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
	@$(RM) ccdftt
	@$(RM) ccdftt.exe
	@$(RM) p
	@$(RM) pp
	@$(RM) p.cpp
