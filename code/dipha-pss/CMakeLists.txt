#  Copyright 2014 IST Austria
#
#  Contributed by: Ulrich Bauer, Jan Reininghaus
#
#  This file is part of DIPHA.
#
#  DIPHA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  DIPHA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with DIPHA.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

project("Dipha-PSS")

set(CMAKE_BUILD_TYPE "Release")

if(MSVC OR MSVC_IDE)
    if( MSVC_VERSION LESS 1700 )
        message(FATAL_ERROR "The project requires C++11 features. " 
        "You need at least Visual Studio 11 (Microsoft Visual Studio 2012)")
    endif()        
    
    # Compilation using Microsoft Visual Studio 2013 seems to be stuck in an infinite loop ...
    if( MSVC_VERSION GREATER 1700 )
        set(CMAKE_GENERATOR_TOOLSET "v110" CACHE STRING "Platform Toolset" FORCE)               
    endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_definitions(-std=c++11 -w)
endif()

# Let user decide whether we want to build FGT support
option(USE_FGT "Use FGT support." OFF)

# Per default, no paths are set
set(FIGTREE_INCLUDE "" CACHE FILEPATH "Path to figtree include directory")
set(FIGTREE_LIB "" CACHE FILEPATH "Path to figtree library")

# In case we want FGT support, we have to set the lib/include dirs
if (USE_FGT)
	if (FIGTREE_LIB STREQUAL "")
		message(FATAL_ERROR "Variable FIGTREE_LIB not set.")
	endif()
	if (FIGTREE_INCLUDE STREQUAL "")
		message(FATAL_ERROR "Variable FIGTREE_INCLUDE not set.")
	endif ()
endif ()


FIND_PACKAGE(MPI REQUIRED)

if (USE_FGT)
	include_directories (../dipha/include ../dipha/externals ${MPI_INCLUDE_PATH} ${FIGTREE_INCLUDE})
	link_directories (${FIGTREE_LIB})
else(USE_FGT)
	include_directories (../dipha/include ../dipha/externals ${MPI_INCLUDE_PATH})
endif(USE_FGT)
	
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

FILE(GLOB_RECURSE all_dipha_includes "../dipha/include/dipha/*.h")
FILE(GLOB dipha_general_includes "../dipha/include/dipha/*.h")
FILE(GLOB dipha_data_structures_includes "../dipha/include/dipha/data_structures/*.h")
FILE(GLOB dipha_input_includes "../dipha/include/dipha/inputs/*.h")
FILE(GLOB dipha_output_includes "../dipha/include/dipha/outputs/*.h")
FILE(GLOB dipha_algorithms_includes "../dipha/include/dipha/algorithms/*.h")

add_executable (diagram_distance src/diagram_distance.cpp)
if (USE_FGT)
	target_link_libraries(diagram_distance ${MPI_LIBRARIES} figtree ann_figtree_version)
	target_compile_definitions(diagram_distance PRIVATE USE_FGT=1)
else (USE_FGT)
	target_link_libraries(diagram_distance ${MPI_LIBRARIES})
endif(USE_FGT)

source_group(Header\ Files FILES ${dipha_general_includes})
source_group(Header\ Files\\inputs FILES ${dipha_input_includes})
source_group(Header\ Files\\outputs FILES ${dipha_output_includes})
source_group(Header\ Files\\data_structures FILES ${dipha_data_structures_includes})
source_group(Header\ Files\\algorithms FILES ${dipha_algorithms_includes})
   
