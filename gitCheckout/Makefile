
# Makefile for G0 Analysis Code
# (modified version of $(ROOTSYS)/test/Makefile)
#
# March 3, 2000  Author: W. Korsch
# April 25, 2000 	 E. Beise, modify to work with Alpha DUX4.0

ObjSuf        = o
SrcSuf        = C
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o  
OFFICIAL      =UCNA_Official_Replay

EXE = analyzer

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)   -lNew -lMinuit -lSpectrum 
ROOTGLIBS     = $(shell root-config --glibs)

MYOS := $(subst -,,$(shell uname))
# Linux with egcs	
ifeq ($(MYOS),Linux)
	CXX           = gcc
	CXXFLAGS      = -O3 -Wall -fPIC -g
	LD            = gcc
	LDFLAGS       = -g
	SOFLAGS       = -shared
        INCLUDEPATH   = -Iincludes -I${OFFICIAL}/IOUtils -I${OFFICIAL}/RootUtils -I${OFFICIAL}/BaseTypes -I${OFFICIAL}/MathUtils -I${OFFICIAL}/Calibration -I${OFFICIAL}/Analysis -I${OFFICIAL}/Studies -I${OFFICIAL}/Physics -I/usr/include/mysql
	CXXFLAGS     += $(ROOTCFLAGS) $(INCLUDEPATH)
	LIBS          = $(ROOTLIBS) -L${OFFICIAL} #-lUCNA
	GLIBS         = $(ROOTGLIBS)
	LIBS	     += -lstdc++ -lz -lRMySQL
endif

#####################################################################
#
#  Generic compilation and linking step to make an executable from
#  a single *.C file
#
%: src/%.$(SrcSuf) %.h
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $< $(LIBS) $(GLIBS)

OBJ     =  src/analysis.o src/Run.o src/Beta_Run.o src/Bck_Run.o src/Octet.o src/analysisMathFunc.o 

clean:
		@rm -f $(OBJ) *~  core

cleanall:
		@rm -f $`(OBJ) $(EXE) *~  core

INCPATH ="includes"

all: analyzer

analyzer: $(OBJ)
	$(CXX) $(CXXFLAGS) -o analyzer $(LIBS) $(OBJ) src/analysisPlottingFunc.C -L $(LIBS) $(GLIBS)	






