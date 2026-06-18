OBJECTS = Job.o JobManager.o Jobs/ComputeCondensedLinearResponse.o Jobs/ComputeElectronDensity.o Jobs/ComputeDescriptors.o Jobs/ComputeEnergyWithPointCharges.o Jobs/ComputeGridDifference.o Jobs/ComputeIntegrals.o Jobs/ComputePartialCharges.o Jobs/ConvertOrbitals.o Jobs/Help.o Jobs/MakeDensityCube.o Jobs/MakeELFCube.o Jobs/MakeOrbitalsCube.o Jobs/RunLambdaDiagnostic.o


include $(LIBCDFTTDIR)/CONFIG

CFLAGS =  $(COMMONCFLAGS) 

all: compiling 

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

Jobs/%.o : Jobs/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

compiling: $(OBJECTS)

clean:
	@$(RM) *o *.out *.log *.inp *~ 
	@$(RM) Jobs/*.o

#make dependance
SRCS = $(OBJECTS:%.o=%.cpp)
dep :
	$(CC)  $(CFLAGS) -MM $(SRCS) > Dep.mk

-include Dep.mk

