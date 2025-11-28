# CMakeLists.txt in main

set (HFILES
    TTMDecay.h
    TTMDecayChannel.h
    TTMDensObj.h
    TTMIDObj.h
    TTMParameter.h
    TTMParameterSet.h
    TTMParameterSetBQ.h
    TTMParameterSetBSQ.h
    TTMParameterSetCanBSQ.h
    TTMParameterSetSQ.h
    TTMParticle.h
    TTMParticleSet.h
    TTMYield.h
    TThermalFit.h
    TThermalFitBQ.h
    TThermalFitBSQ.h
    TThermalFitCanBSQ.h
    TThermalFitSQ.h
    TThermalModel.h
    TThermalModelBQ.h
    TThermalModelBSQ.h
    TThermalModelCanBSQ.h
    TThermalModelSQ.h
    TThermalParticle.h
    TThermalParticleBQ.h
    TThermalParticleBSQ.h
    TThermalParticleCanBSQ.h
    TThermalParticleSQ.h
)

set(SRCFILES
    TTMDecay.cxx
    TTMDecayChannel.cxx
    TTMDensObj.cxx
    TTMIDObj.cxx
    TTMParameter.cxx
    TTMParameterSet.cxx
    TTMParameterSetBQ.cxx
    TTMParameterSetBSQ.cxx
    TTMParameterSetCanBSQ.cxx
    TTMParameterSetSQ.cxx
    TTMParticle.cxx
    TTMParticleSet.cxx
    TTMYield.cxx
    TThermalFit.cxx
    TThermalFitBQ.cxx
    TThermalFitBSQ.cxx
    TThermalFitCanBSQ.cxx
    TThermalFitSQ.cxx
    TThermalModel.cxx
    TThermalModelBQ.cxx
    TThermalModelBSQ.cxx
    TThermalModelCanBSQ.cxx
    TThermalModelSQ.cxx
    TThermalParticle.cxx
    TThermalParticleBQ.cxx
    TThermalParticleBSQ.cxx
    TThermalParticleCanBSQ.cxx
    TThermalParticleSQ.cxx
)

ROOT_STANDARD_LIBRARY_PACKAGE(
        THERMUS
    HEADERS
        ${HFILES}
    SOURCES
        ${SRCFILES}
    DEPENDENCIES
        ${ROOT_LIBRARIES})

set_target_properties(THERMUS PROPERTIES PUBLIC_HEADER "${HFILES}")

install(TARGETS THERMUS
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/main/libTHERMUS_rdict.pcm 
    ${CMAKE_BINARY_DIR}/main/libTHERMUS.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
