# CMakeLists.txt in main

FILE(GLOB SRCFILES *.cxx)
FILE(GLOB HFILES TT*.h)

ROOT_GENERATE_DICTIONARY(ThermusDict ${HFILES} LINKDEF THERMUSClassesLinkDef.h)

add_library(THERMUS SHARED
 # ${SRCFILES}
 ../functions/BQConstrainBDens.cxx
 ../functions/BQConstrainQEN.cxx
 ../functions/BQConstrainST3.cxx
 ../functions/BQfuncQBDens.cxx
 ../functions/BQfuncQST3.cxx
 ../functions/BQConstrainEN.cxx
 ../functions/BQConstrainQNetBDens.cxx
 ../functions/BQfuncBDens.cxx
 ../functions/BQfuncQEN.cxx
 ../functions/BQfuncST3.cxx
 ../functions/BQConstrainQ.cxx
 ../functions/BQConstrainQPercolation.cxx
 ../functions/BQfuncEN.cxx
 ../functions/BQfuncQNetBDens.cxx
 ../functions/BQConstrainQBDens.cxx
 ../functions/BQConstrainQST3.cxx
 ../functions/BQfuncQ.cxx
 ../functions/BQfuncQPercolation.cxx
 ../functions/brent.cxx
 ../functions/broyden.cxx
 ../functions/broydens.cxx
 ../functions/BSQCanonicalFncs.cxx
 ../functions/BSQConstrainSQBDens.cxx
 ../functions/BSQfuncQ.cxx
 ../functions/BSQfuncSQCb.cxx
 ../functions/BSQConstrainBSQ.cxx
 ../functions/BSQConstrainSQC.cxx
 ../functions/BSQfuncQQ.cxx
 ../functions/BSQfuncSQEN.cxx
 ../functions/BSQConstrainQ.cxx
 ../functions/BSQConstrainSQCb.cxx
 ../functions/BSQfuncS.cxx
 ../functions/BSQfuncSQNetBDens.cxx
 ../functions/BSQConstrainQQ.cxx
 ../functions/BSQConstrainSQEN.cxx
 ../functions/BSQfuncSBDens.cxx
 ../functions/BSQfuncSQPercolation.cxx
 ../functions/BSQConstrainS.cxx
 ../functions/BSQConstrainSQNetBDens.cxx
 ../functions/BSQfuncSEN.cxx
 ../functions/BSQfuncSQST3.cxx
 ../functions/BSQConstrainSBDens.cxx
 ../functions/BSQConstrainSQPercolation.cxx
 ../functions/BSQfuncSPercolation.cxx
 ../functions/BSQfuncSST3.cxx
 ../functions/BSQConstrainSEN.cxx
 ../functions/BSQConstrainSQST3.cxx
 ../functions/BSQfuncSQ.cxx
 ../functions/BSQConstrainSPercolation.cxx
 ../functions/BSQConstrainSST3.cxx
 ../functions/BSQfuncSQBDens.cxx
 ../functions/BSQConstrainSQ.cxx
 ../functions/BSQfuncBSQ.cxx
 ../functions/BSQfuncSQC.cxx
 ../functions/ExclVolPressureFunc.cxx
 ../functions/FindExclVolPressure.cxx
 ../functions/fit_function.cxx
 ../functions/Functions.cxx
 ../functions/Int_2_String.cxx
 ../functions/IntegrationRoutines.cxx
 ../functions/IntegrationRoutines2D.cxx
 ../functions/Minuit_function.cxx
 ../main/TTMDecay.cxx
 ../main/TTMDecayChannel.cxx
 ../main/TTMDensObj.cxx
 ../main/TTMIDObj.cxx
 ../main/TTMParameter.cxx
 ../main/TTMParameterSet.cxx
 ../main/TTMParameterSetBQ.cxx
 ../main/TTMParameterSetBSQ.cxx
 ../main/TTMParameterSetCanBSQ.cxx
 ../main/TTMParameterSetSQ.cxx
 ../main/TTMParticle.cxx
 ../main/TTMParticleSet.cxx
 ../main/TTMYield.cxx
 ../main/TThermalFit.cxx
 ../main/TThermalFitBQ.cxx
 ../main/TThermalFitBSQ.cxx
 ../main/TThermalFitCanBSQ.cxx
 ../main/TThermalFitSQ.cxx
 ../main/TThermalModel.cxx
 ../main/TThermalModelBQ.cxx
 ../main/TThermalModelBSQ.cxx
 ../main/TThermalModelCanBSQ.cxx
 ../main/TThermalModelSQ.cxx
 ../main/TThermalParticle.cxx
 ../main/TThermalParticleBQ.cxx
 ../main/TThermalParticleBSQ.cxx
 ../main/TThermalParticleCanBSQ.cxx
 ../main/TThermalParticleSQ.cxx
 ../main/ThermusDict.cxx
 )

set_target_properties(THERMUS PROPERTIES
VERSION ${PROJECT_VERSION}
)
target_link_libraries(THERMUS ${ROOT_LIBRARIES} ROOT::Minuit GSL::gsl GSL::gslcblas)

