# CMakeLists.txt in main

FILE(GLOB SRCFILES *.cxx)
FILE(GLOB HFILES TT*.h)

#add_library(mylib some.c another.c)
#set_target_properties(mylib PROPERTIES PUBLIC_HEADER "some.h;another.h")
#INSTALL(TARGETS mylib
#        LIBRARY DESTINATION some/libpath
#        PUBLIC_HEADER DESTINATION some/includepath
#)

ROOT_GENERATE_DICTIONARY(ThermusDict ${HFILES} LINKDEF THERMUSClassesLinkDef.h)

add_library(THERMUS SHARED ${SRCFILES} ThermusDict.cxx)
set_target_properties(THERMUS PROPERTIES PUBLIC_HEADER "${HFILES}")
target_link_libraries(THERMUS ${ROOT_LIBRARIES} GSL::gsl GSL::gslcblas)

install(TARGETS THERMUS)
#install(TARGETS )
