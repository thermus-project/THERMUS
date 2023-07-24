#########################
ARCH          = $(shell root-config --arch)

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


CXX           = g++
CXXFLAGS      = -O  -fPIC
LD            = g++
LDFLAGS       = -O  
ifeq ($(ARCH),macosx64)
SOFLAGS       = -m64 -dynamiclib -single_module -undefined dynamic_lookup
else
SOFLAGS       = -shared
endif


CXXFLAGS     += $(ROOTCFLAGS)
CXXFLAGS     += -I$(THERMUS)/main -I$(THERMUS)/include -I$(THERMUS)/nrc
LIBS          = $(ROOTLIBS) $(SYSLIBS) -L/usr/local/lib -lgsl -lgslcblas -lm
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------


SEP_CLASSESH	= $(THERMUS)/main/TTMParticle.h \
		  $(THERMUS)/main/TTMDecay.h \
		  $(THERMUS)/main/TTMDecayChannel.h \
		  $(THERMUS)/main/TTMParticleSet.h \
		  $(THERMUS)/main/TTMParameter.h \
		  $(THERMUS)/main/TTMParameterSet.h \
		  $(THERMUS)/main/TTMParameterSetBSQ.h \
		  $(THERMUS)/main/TTMParameterSetBQ.h \
		  $(THERMUS)/main/TTMParameterSetCanBSQ.h \
		  $(THERMUS)/main/TThermalParticle.h \
		  $(THERMUS)/main/TThermalParticleBSQ.h \
		  $(THERMUS)/main/TThermalParticleBQ.h \
		  $(THERMUS)/main/TThermalParticleCanBSQ.h \
		  $(THERMUS)/main/TThermalModel.h \
		  $(THERMUS)/main/TTMYield.h \
		  $(THERMUS)/main/TTMDensObj.h \
		  $(THERMUS)/main/TThermalModelBSQ.h \
		  $(THERMUS)/main/TThermalModelBQ.h \
		  $(THERMUS)/main/TThermalModelCanBSQ.h \
		  $(THERMUS)/main/TThermalFit.h \
		  $(THERMUS)/main/TThermalFitBSQ.h \
		  $(THERMUS)/main/TThermalFitBQ.h \
		  $(THERMUS)/main/TThermalFitCanBSQ.h \
		  $(THERMUS)/main/TTMIDObj.h \

FNCSO        	= BQConstrainQ.o \
		  BSQConstrainS.o \
		  BSQConstrainQ.o \
		  BSQConstrainSQ.o \
		  BSQConstrainSQC.o \
		  BSQConstrainSQCb.o \
		  BSQConstrainBSQ.o \
		  BSQConstrainQQ.o \
		  BSQConstrainSQEN.o \
		  BSQConstrainSEN.o \
		  BSQConstrainSQST3.o \
		  BSQConstrainSST3.o \
		  BSQConstrainSQPercolation.o \
		  BSQConstrainSPercolation.o \
		  BQConstrainEN.o \
		  BQConstrainQEN.o \
		  BQConstrainST3.o \
		  BQConstrainBDens.o \
		  BQConstrainQST3.o \
	  	  BQConstrainQBDens.o \
	  	  BQConstrainQNetBDens.o \
		  BQConstrainQPercolation.o \
		  BSQConstrainSBDens.o \
		  BSQConstrainSQBDens.o \
		  BSQConstrainSQNetBDens.o \
		  BQfuncEN.o \
          BQfuncQEN.o \
          BQfuncST3.o \
          BQfuncBDens.o \
          BQfuncQST3.o \
          BQfuncQBDens.o \
          BQfuncQNetBDens.o \
		  BQfuncQ.o \
		  BQfuncQPercolation.o \
		  BSQfuncSQPercolation.o \
		  BSQfuncSPercolation.o \
	      BSQfuncS.o \
          BSQfuncQ.o \
		  BSQfuncSQ.o \
		  BSQfuncSQC.o \
		  BSQfuncSQCb.o \
		  BSQfuncBSQ.o \
		  BSQfuncQQ.o \
		  BSQfuncSQEN.o \
		  BSQfuncSEN.o \
		  BSQfuncSQST3.o \
		  BSQfuncSST3.o \
		  BSQfuncSBDens.o \
		  BSQfuncSQBDens.o \
		  BSQfuncSQNetBDens.o \
		  fit_function.o \
		  Minuit_function.o \
		  Functions.o \
		  BSQCanonicalFncs.o \
		  IntegrationRoutines.o \
		  IntegrationRoutines2D.o \
		  Int_2_String.o \
		  FindExclVolPressure.o \
          ExclVolPressureFunc.o \
          broyden.o \
          broydens.o \
          brent.o


FNCSS   = BQConstrainQ.cxx \
		  BSQConstrainS.cxx \
		  BSQConstrainQ.cxx \
		  BSQConstrainSQ.cxx \
		  BSQConstrainSQC.cxx \
		  BSQConstrainSQCb.cxx \
		  BSQConstrainBSQ.cxx \
		  BSQConstrainQQ.cxx \
		  BSQConstrainSQEN.cxx \
		  BSQConstrainSEN.cxx \
		  BSQConstrainSQST3.cxx \
		  BSQConstrainSST3.cxx \
		  BSQConstrainSQPercolation.cxx \
		  BSQConstrainSPercolation.cxx \
		  BQConstrainEN.cxx \
          BQConstrainQEN.cxx \
          BQConstrainST3.cxx \
          BQConstrainBDens.cxx \
          BQConstrainQST3.cxx \
          BQConstrainQBDens.cxx \
          BQConstrainQNetBDens.cxx \
		  BQConstrainQPercolation.cxx \
		  BSQConstrainSBDens.cxx \
		  BSQConstrainSQBDens.cxx \
		  BSQConstrainSQNetBDens.cxx \
		  BQfuncEN.cxx \
		  BQfuncQEN.cxx \
		  BQfuncST3.cxx \
		  BQfuncBDens.cxx \
		  BQfuncQST3.cxx \
		  BQfuncQBDens.cxx \
		  BQfuncQNetBDens.cxx \
		  BQfuncQ.cxx \
		  BQfuncQPercolation.cxx \
	      BSQfuncS.cxx \
	      BSQfuncQ.cxx \
		  BSQfuncSQ.cxx \
		  BSQfuncSQC.cxx \
		  BSQfuncSQCb.cxx \
		  BSQfuncBSQ.cxx \
		  BSQfuncQQ.cxx \
		  BSQfuncSQEN.cxx \
		  BSQfuncSEN.cxx \
		  BSQfuncSQST3.cxx \
		  BSQfuncSST3.cxx \
		  BSQfuncSBDens.cxx \
		  BSQfuncSQBDens.cxx \
		  BSQfuncSQNetBDens.cxx \
		  BSQfuncSQPercolation.cxx \
		  BSQfuncSPercolation.cxx \
		  fit_function.cxx \
		  Minuit_function.cxx \
		  Functions.cxx \
		  BSQCanonicalFncs.cxx \
		  IntegrationRoutines.cxx \
		  IntegrationRoutines2D.cxx \
		  Int_2_String.cxx \
		  FindExclVolPressure.cxx \
		  ExclVolPressureFunc.cxx \
          broyden.cxx \
          broydens.cxx \
          brent.cxx

FNCSSO		= $(THERMUS)/lib/libFunctions.$(DllSuf)

OBJS          	= $(FNCSO) 
PROGRAMS      	= $(FNCSSO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   Functions

all:            $(PROGRAMS)

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

###

Functions:	$(FNCSSO)	

$(FNCSSO):	$(FNCSO)
	$(LD) $(SOFLAGS) $(LIBS) $(LDFLAGS) $^ $(OutPutOpt) $@


$(FNCSO): 	$(FNCSS) $(SEP_CLASSESH) 
	$(CXX) $(CXXFLAGS) -c $(FNCSS) 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
