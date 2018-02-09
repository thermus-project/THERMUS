# CMakeLists.txt in main

FILE(GLOB SRCFILES *.cxx)
FILE(GLOB HFILES TT*.h)

ROOT_GENERATE_DICTIONARY(ThermusDict ${HFILES} LINKDEF THERMUSClassesLinkDef.h)

add_library(THERMUS SHARED ${SRCFILES} ThermusDict.cxx)
target_link_libraries(THERMUS ${ROOT_LIBRARIES} GSL::gsl GSL::gslcblas)
