cmake_minimum_required (VERSION 3.21 FATAL_ERROR)
project (Thermus VERSION 4.0.0)

# Find GSL
FIND_PACKAGE( GSL REQUIRED )
if ( NOT GSL_FOUND )
message(FATAL_ERROR "GSL is required but missing !")
endif ( NOT GSL_FOUND )

# Find xxHash FIXME and liburing2
#FIND_PACKAGE( xxHash REQUIRED )
#if ( NOT xxHash_FOUND )
#message(FATAL_ERROR "xxHash is required but missing !")
#endif ( NOT xxHash_FOUND )

# Check ROOT and Mathmore
#FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathMore )
FIND_PACKAGE( ROOT CONFIG REQUIRED COMPONENTS MathMore )
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )

include(${ROOT_USE_FILE})

#GNU like structure
include(GNUInstallDirs)
# Put libs under Thermus in libdir
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR}/Thermus)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/Thermus)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

include_directories("include" "main" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

if (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)
set(CMAKE_SHARED_LIBRARY_SUFFIX .so)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)

if(${CMAKE_SYSTEM} MATCHES Darwin)
  set (CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup -lm")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

message(STATUS "Customised the linking flags depending on the platform: ${CMAKE_SHARED_LINKER_FLAGS}")

foreach(p
 CMP0042 # MACOSX_RPATH is enabled by default (v3.0)
 CMP0053 # Simplify variable reference and escape sequence evaluation (v3.1.3)
)
if ( POLICY ${p} )
cmake_policy(SET ${p} NEW)
endif()
endforeach()

message(STATUS "Source folder: ${Thermus_SOURCE_DIR}")
message(STATUS "               which should correspond to $THERMUS ")
message(STATUS "CMake platform: ${CMAKE_SYSTEM}")
message(STATUS "Build folder: ${Thermus_BINARY_DIR}")
#message(STATUS "Shared Object Library folder: ${Thermus_BINARY_DIR}/lib")

add_subdirectory(functions)
add_subdirectory(main)
add_subdirectory(particles)
# Don't include do if WITHOUT_TEX_DOC is set
if(NOT WITHOUT_TEX_DOC)
  add_subdirectory(doc)
endif()
add_subdirectory(scripts)
add_subdirectory(tests)

# Also publish the README file
set(README_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)

install(FILES ${README_FILE} DESTINATION ${CMAKE_INSTALL_DOCDIR})
