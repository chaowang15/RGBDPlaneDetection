VERSION = MRF2.2

SRC =  mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp \
       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp \
       TRW-S.cpp BP-S.cpp

CC = g++

WARN = -W -Wall
OPT ?= -O3
CPPFLAGS = $(OPT) $(WARN) -DUSE_64_BIT_PTR_CAST
#CPPFLAGS = $(OPT) $(WARN)   ### use this line instead to compile on 32-bit systems

OBJ = $(SRC:.cpp=.o)

all: libMRF.a example

libMRF.a: $(OBJ)
	rm -f libMRF.a
	ar ruc libMRF.a $(OBJ)
	ranlib libMRF.a

example: libMRF.a example.cpp
	$(CC) -o example example.cpp -L. -lMRF

clean: 
	rm -f $(OBJ) core core.* *.stackdump *.bak

allclean: clean
	rm -f libMRF.a example example.exe

depend:
	@makedepend -Y -- $(CPPFLAGS) -- $(SRC) 2>> /dev/null

# DO NOT DELETE THIS LINE -- make depend depends on it.

mrf.o: mrf.h
ICM.o: ICM.h mrf.h LinkedBlockList.h
GCoptimization.o: energy.h graph.h block.h mrf.h GCoptimization.h
GCoptimization.o: LinkedBlockList.h
graph.o: graph.h block.h mrf.h
maxflow.o: graph.h block.h mrf.h
MaxProdBP.o: MaxProdBP.h mrf.h LinkedBlockList.h regions-new.h
LinkedBlockList.o: LinkedBlockList.h
regions-maxprod.o: MaxProdBP.h mrf.h LinkedBlockList.h regions-new.h
TRW-S.o: TRW-S.h mrf.h typeTruncatedQuadratic2D.h
BP-S.o: BP-S.h mrf.h typeTruncatedQuadratic2D.h
