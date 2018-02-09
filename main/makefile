######################
ARCH          = $(shell root-config --arch)

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
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
CXXFLAGS     += -I$(THERMUS)/main -I$(THERMUS)/include
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------


SEP_CLASSESO    = $(THERMUS)/main/TTMParticle.$(ObjSuf) \
		  $(THERMUS)/main/TTMDecay.$(ObjSuf) \
		  $(THERMUS)/main/TTMParticleSet.$(ObjSuf) \
		  $(THERMUS)/main/TTMParameter.$(ObjSuf) \
		  $(THERMUS)/main/TTMParameterSet.$(ObjSuf) \
		  $(THERMUS)/main/TTMParameterSetBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TTMParameterSetBQ.$(ObjSuf) \
		  $(THERMUS)/main/TTMParameterSetCanBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalParticle.$(ObjSuf) \
		  $(THERMUS)/main/TThermalParticleBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalParticleBQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalParticleCanBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalModel.$(ObjSuf) \
		  $(THERMUS)/main/TTMYield.$(ObjSuf) \
		  $(THERMUS)/main/TTMDensObj.$(ObjSuf) \
		  $(THERMUS)/main/TThermalModelBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalModelBQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalModelCanBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalFit.$(ObjSuf) \
		  $(THERMUS)/main/TThermalFitBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalFitBQ.$(ObjSuf) \
		  $(THERMUS)/main/TThermalFitCanBSQ.$(ObjSuf) \
		  $(THERMUS)/main/TTMIDObj.$(ObjSuf) \
		  $(THERMUS)/main/TTMDecayChannel.$(ObjSuf) \

CLASSDICTO	= $(THERMUS)/main/ClassesDict.$(ObjSuf)

SEP_CLASSESS    = $(THERMUS)/main/TTMParticle.$(SrcSuf) \
		  $(THERMUS)/main/TTMDecay.$(SrcSuf) \
		  $(THERMUS)/main/TTMParticleSet.$(SrcSuf) \
		  $(THERMUS)/main/TTMParameter.$(SrcSuf) \
		  $(THERMUS)/main/TTMParameterSet.$(SrcSuf) \
		  $(THERMUS)/main/TTMParameterSetBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TTMParameterSetBQ.$(SrcSuf) \
		  $(THERMUS)/main/TTMParameterSetCanBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalParticle.$(SrcSuf) \
		  $(THERMUS)/main/TThermalParticleBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalParticleBQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalParticleCanBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalModel.$(SrcSuf) \
		  $(THERMUS)/main/TTMYield.$(SrcSuf) \
		  $(THERMUS)/main/TTMDensObj.$(SrcSuf) \
		  $(THERMUS)/main/TThermalModelBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalModelBQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalModelCanBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalFit.$(SrcSuf) \
		  $(THERMUS)/main/TThermalFitBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalFitBQ.$(SrcSuf) \
		  $(THERMUS)/main/TThermalFitCanBSQ.$(SrcSuf) \
		  $(THERMUS)/main/TTMIDObj.$(SrcSuf) \
		  $(THERMUS)/main/TTMDecayChannel.$(SrcSuf) \

CLASSDICTS	= $(THERMUS)/main/ClassesDict.$(SrcSuf)

SEP_CLASSESH	= $(THERMUS)/main/TTMParticle.h \
		  $(THERMUS)/main/TTMDecay.h \
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
		  $(THERMUS)/main/TTMDecayChannel.h \

THERMALSO	= $(THERMUS)/lib/libTHERMUS.$(DllSuf)

OBJS          	= $(SEP_CLASSESO) $(CLASSDICTO)  
PROGRAMS      	= $(THERMALSO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   Classes Thermal

all:            $(PROGRAMS)
		@echo "Copying dictionary $@..."
		@cp -a -f *.pcm $(THERMUS)/lib

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository
		@rm -f $(THERMUS)/lib/*.pcm

###

Thermal:	$(THERMALSO)	

$(THERMALSO):	$(SEP_CLASSESO) $(CLASSDICTO) 
	$(LD) $(SOFLAGS) $(LIBS) $(LDFLAGS) $^ $(OutPutOpt) $@

Classes:	$(SEP_CLASSESO) $(CLASSDICTO)

$(SEP_CLASSESO): 	$(SEP_CLASSESS) $(SEP_CLASSESH) 
	$(CXX) $(CXXFLAGS) -c $(SEP_CLASSESS) 

$(CLASSDICTO): $(CLASSDICTS)
	$(CXX) $(CXXFLAGS) -c $(CLASSDICTS) 

$(CLASSDICTS): $(SEP_CLASSESH) $(THERMUS)/main/THERMUSClassesLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcling -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
