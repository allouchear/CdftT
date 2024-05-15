include CONFIG


ANALYTIC = analytic
NUMERIC = numeric
COMMON  = common

ANALYTICOBJ = $(ANALYTIC)/*.o
NUMERICOBJ  = $(NUMERIC)/*.o
COMMONSOBJ  = $(COMMON)/*.o

SUBDIRS = \
	$(COMMON) $(NUMERIC) $(ANALYTIC)

OBJECTS = \
	$(ANALYTICOBJ) $(NUMERICOBJ) $(COMMONSOBJ)

TARGETS = all clean dep
 
.PHONY: $(SUBDIRS)

$(TARGETS): $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --directory=$@ $(MAKECMDGOALS)

all: cchemilib

cchemilib: $(SUBDIRS)
	$(CC) -fPIC -shared $(LDFLAGS) $(OBJECTS) -o ../../lib/libcdftt.so  $(WIN32LIB)  $(OMPLIB) $(CLLIB) -g -lm
 
clean:
	@$(RM) p
	@$(RM) pp
	@$(RM) p.cpp
	@$(RM) *.log
	@$(RM) *.out
	@$(RM) *.txt
