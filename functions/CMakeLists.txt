# CMakeLists.txt in functions

FILE(GLOB SRCFILES *.cxx)

add_library(Functions SHARED ${SRCFILES})
target_link_libraries(Functions ROOT::Minuit GSL::gsl GSL::gslcblas)

install(TARGETS Functions)
