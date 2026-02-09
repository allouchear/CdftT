OBJECTS = ExcitedState.o Orbitals.o OrbitalsSave.o SlaterDeterminant.o


include $(LIBCDFTTDIR)/CONFIG

CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

profile: CFLAGS = $(COMMONCFLAGS) $(PROFILEFLAGS)
profile: compiling

clean:
	@$(RM) *o *.out *.log *.inp *~ *.mk

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

