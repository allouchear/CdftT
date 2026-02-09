OBJECTS = GridPoints6.o GridPoints14.o GridPoints26.o GridPoints38.o GridPoints50.o GridPoints74.o GridPoints86.o GridPoints110.o GridPoints146.o GridPoints170.o GridPoints194.o GridPoints230.o GridPoints266.o GridPoints302.o GridPoints350.o GridPoints434.o GridPoints590.o GridPoints770.o GridPoints974.o GridPoints1202.o GridPoints1454.o GridPoints1730.o GridPoints2030.o GridPoints2354.o GridPoints5810.o GridPoints.o Becke.o


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

