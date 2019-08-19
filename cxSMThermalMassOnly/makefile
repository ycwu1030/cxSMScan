
# LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
# LT = /Users/mac/work/LoopTools/x86_64-Darwin/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
# CXX = $(LT)/../bin/f++
CXX = g++
# FFLAG = -I$(LT)/../include -I$(INCDIR) #-DDEBUG
# FLIBS = -L$(LT) -looptools

gsllibs = $(shell gsl-config --libs)
gslcflags = $(shell gsl-config --cflags)
#boostlibs = -L/usr/local/Cellar/boost/1.70.0/lib -lboost_math_c99l
#boostcflags = -I/usr/local/Cellar/boost/1.70.0/include

FFLAG = -I$(INCDIR) $(gslcflags) $(boostcflags)
FLIBS = $(gsllibs) $(boostlibs)

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

all: cxSMCritical_Scan.x

.PHONY: clean

.SECONDARY: $(OBJ)

%.x:%.cpp $(OBJ)
	$(CXX) $(FFLAG) -o $@ $< $(OBJ) $(FLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(FFLAG) -c $< -o $@


clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

