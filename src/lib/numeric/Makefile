OBJECTS = Domain.o Grid.o functions.o GridCP.o


include $(LIBCDFTTDIR)/CONFIG

CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

