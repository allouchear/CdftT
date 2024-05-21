include ../../../CONFIG

CFLAGS =  $(COMMONCFLAGS) 

ANALYTIC = ../../lib/analytic
NUMERIC = ../../lib/numeric
COMMON  = ../../lib/common

ANALYTICOBJ = $(ANALYTIC)/*.o
NUMERICOBJ  = $(NUMERIC)/*.o
COMMONSOBJ  = $(COMMON)/*.o

OBJECTS = main.o \
	$(ANALYTICOBJ) $(NUMERICOBJ) $(COMMONSOBJ)

all: compiling  cdftt

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) cdftt *o *.out *.log *.inp *~ 

cdftt: 
	$(CC) $(LDFLAGS) $(OBJECTS) -o cdftt $(OMPLIB) -g -lm -lgomp 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk



-include Dep.mk

