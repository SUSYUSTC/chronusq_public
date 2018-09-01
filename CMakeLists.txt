#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2018 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#
cmake_minimum_required(VERSION 3.0)
project(ChronusQ CXX C Fortran)
set(ChronusQ_VERSION_MAJOR 0)
set(ChronusQ_VERSION_MINOR 2)
set(ChronusQ_VERSION_PATCH 0)


# Options
option(CQ_ENABLE_MPI       "Enable MPI parallelism"          OFF)
#option(CQ_ENABLE_SCALAPACK "Enable ScaLAPACK"                OFF)
#option(CQ_ENABLE_CXXBLACS  "Enable CXXBLACS"                 OFF)
option(ENABLE_COVERAGE     "Enable coverage and profiling"   OFF)
option(FORCE_SHARED_BOOST  "Force link to shared Boost Libs" OFF)
option(CQ_LINALG_USESYSTEM "Use System defaults for LA libs" OFF)
option(FORCE_SHARED_BOOST  "Force link to shared Boost Libs" OFF)


# Add specified paths to CMAKE_PREFIX_PATH
if( MPI_PREFIX )
  set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MPI_PREFIX} )
endif()

# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Intel C++11
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()



# Parallelism

message( "\n\n" )
message( "ChronusQ Parallelism Settings" )

# OpenMP
message( "" )
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if( MPI_PREFIX )
  # MPI_PREFIX paths to CMAKE_PREFIX_PATH
  message( "" )
  message( STATUS "Adding ${MPI_PREFIX} to CMAKE_PREFIX_PATH" )
  set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MPI_PREFIX} )
endif()

# MPI
if(CQ_ENABLE_MPI)
  include(HandleMPI)
endif()


message( "\n\n" )









# FLAGS
include(HandleCXXFlags)


# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

# Linear Algebra
include(HandleLA)



# Non-trivial dependencies
message( "ChronusQ Dependencies" )


include(FindLibint)
include(FindLibXC)
include(Boost)
include(pcm)
include(cnpy)

set(CQEX_LINK ${CQ_EXT_LINK} )

# HDF5
find_package(HDF5 REQUIRED CXX)
include_directories("${HDF5_INCLUDE_DIRS}")
list(APPEND CQEX_LINK ${HDF5_CXX_LIBRARIES} )

message( "\n\n" )




# Compiling with coverage report
if(ENABLE_COVERAGE)

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "ENABLE_COVERAGE Only Compatible with GNU Compilers!")
  else()
    message(STATUS "Enabling Coverage and Profiling -- Warning! May impede performance")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  endif()
endif()



set(BASIS_PATH     "${PROJECT_SOURCE_DIR}/basis"     )
set(AUTHOR_PATH    "${PROJECT_SOURCE_DIR}/AUTHORS"   )

configure_file(
  "${PROJECT_SOURCE_DIR}/include/chronusq_config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/chronusq_config.hpp"
)

# Recurse into subdirectories
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
