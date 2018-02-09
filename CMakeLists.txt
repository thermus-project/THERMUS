cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
project (Thermus)
# checked on 30/01/2018 with root-v6-12-04 (B.H.)

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

set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)
include(RootInstallDirs)

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("include" "main" ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})

if (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)
set(CMAKE_SHARED_LIBRARY_SUFFIX .so)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL Windows)

if(${CMAKE_SYSTEM} MATCHES Darwin)
  set (CMAKE_SHARED_LINKER_FLAGS "-single_module -undefined dynamic_lookup -lm")
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

add_subdirectory(functions)
add_subdirectory(main)
