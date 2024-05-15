main.o: main.cpp
CGTF.o: ../../lib/analytic/CGTF.cpp
GTF.o: ../../lib/analytic/GTF.cpp
MathFunction.o: ../../lib/analytic/MathFunction.cpp
WFX.o: ../../lib/analytic/WFX.cpp
test.o: ../../lib/analytic/test.cpp
Domain.o: ../../lib/numeric/Domain.cpp
Grid.o: ../../lib/numeric/Grid.cpp
test.o: ../../lib/numeric/test.cpp
Atom.o: ../../lib/common/Atom.cpp ../../lib/common/Atom.h \
 ../../lib/common/PeriodicTable.h ../../lib/common/Element.h \
 ../../lib/common/Constants.h
Element.o: ../../lib/common/Element.cpp ../../lib/common/Element.h
PeriodicTable.o: ../../lib/common/PeriodicTable.cpp \
 ../../lib/common/PeriodicTable.h ../../lib/common/Element.h
Structure.o: ../../lib/common/Structure.cpp ../../lib/common/Structure.h \
 ../../lib/common/Atom.h ../../lib/common/PeriodicTable.h \
 ../../lib/common/Element.h
