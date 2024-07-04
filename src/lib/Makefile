include $(LIBCDFTTDIR)/CONFIG


ANALYTIC = analytic
NUMERIC = numeric
COMMON  = common
JOBCONTROL  = jobcontrol

ANALYTICOBJ = $(ANALYTIC)/*/*.o
NUMERICOBJ  = $(NUMERIC)/*.o
COMMONSOBJ  = $(COMMON)/*.o
JOBCONTROLOBJ  = $(JOBCONTROL)/*.o

SUBDIRS = \
	$(COMMON) $(NUMERIC) $(ANALYTIC) $(JOBCONTROL)

OBJECTS = \
	$(ANALYTICOBJ) $(NUMERICOBJ) $(COMMONSOBJ) $(JOBCONTROLOBJ)

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
