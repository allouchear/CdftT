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

all: compiling  

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	gcc  $(CFLAGS) -MM $(SRCS) > Dep.mk

cdftt: 
	$(CC) $(LDFLAGS) $(OBJECTS) -o cdftt $(OMPLIB) -g -lm -lgomp 


-include Dep.mk

