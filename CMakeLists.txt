cmake_minimum_required (VERSION 3.21 FATAL_ERROR)
project (Thermus VERSION 4.0 DESCRIPTION "Thermus v4.0")
# checked on 24/07/2023 with root-v6-28-04 (B.H.)

# Find GSL
FIND_PACKAGE( GSL REQUIRED )
if ( NOT GSL_FOUND )
message(FATAL_ERROR "GSL is required but missing !")
endif ( NOT GSL_FOUND )

# Check ROOT and Mathmore
FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathMore )
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET (CMAKE_INSTALL_BINDIR ${CMAKE_BINARY_DIR}/bin)
SET (CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

include_directories("include" "main" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

if (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)
set(CMAKE_SHARED_LIBRARY_SUFFIX .so)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)

if(${CMAKE_SYSTEM} MATCHES Darwin)
# building a dynamic library containing only one module: single_module
# specifying dynamic lookup for undefined symbols as private definitions:
#
set (CMAKE_SHARED_LINKER_FLAGS "-single_module -undefined dynamic_lookup -lm")
#set (CMAKE_SHARED_LINKER_FLAGS "-single_module -undefined dynamic_lookup --verbose -lm")
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
message(STATUS "Shared Object Library folder: ${Thermus_BINARY_DIR}/lib")

#add_subdirectory(functions)
add_subdirectory(main)
