OBJECTS = Atom.o Element.o PeriodicTable.o Structure.o Descriptors.o


include ../../../CONFIG
CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CC) -g -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

