cmake_minimum_required(VERSION 3.16)

set(CMAKE_Fortran_COMPILER $ENV{HDF5_ROOT}/bin/h5pfc)

project(tristan-ff VERSION 1.0
        DESCRIPTION "Simple relativistic force-free code"
        LANGUAGES Fortran)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(HDF5 REQUIRED COMPONENTS Fortran)

set(DEBUG "OFF" CACHE BOOL "Debug mode")

if (${DEBUG} STREQUAL "ON")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDEBUG")
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fimplicit-none -fbacktrace")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")
  endif()
else()
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -DSoA -xHost -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto")
  endif()
endif()

add_subdirectory(src)
