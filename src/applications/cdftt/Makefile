include $(LIBCDFTTDIR)/CONFIG

CFLAGS =  $(COMMONCFLAGS) 

ANALYTIC = $(LIBCDFTTDIR)/src/lib/analytic
NUMERIC = $(LIBCDFTTDIR)/src/lib/numeric
COMMON  = $(LIBCDFTTDIR)/src/lib/common
JOBCONTROL  = $(LIBCDFTTDIR)/src/lib/jobcontrol

ANALYTICOBJ = $(ANALYTIC)/*/*.o
NUMERICOBJ  = $(NUMERIC)/*.o
COMMONSOBJ  = $(COMMON)/*.o
JOBCONTROLOBJ  = $(JOBCONTROL)/*.o


OBJECTS = cdftt.o \
	$(ANALYTICOBJ) $(NUMERICOBJ) $(COMMONSOBJ) $(JOBCONTROLOBJ)

all: compiling  cdftt

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.inp *~ cdftt

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

cdftt: 
	$(CC) $(LDFLAGS) $(OBJECTS) -o cdftt $(OMPLIB) -g -lm $(ACCLIB)  


-include Dep.mk

