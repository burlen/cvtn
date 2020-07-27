cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(cvt)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Type of build: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if (NOT MSVC)
  if (NOT CMAKE_CXX_FLAGS)
  set(tmp "-fPIC -std=c++17 -Wall -Wextra")
  if (BUILD_STATIC_EXECS)
    set(tmp "${tmp} -static -static-libgcc -static-libstdc++ -pthread -Wl,-Bstatic")
  endif()
  if ((APPLE) AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
    set(tmp "${tmp} -stdlib=libc++")
  endif()
  if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    set(tmp "${tmp} -O3 -march=native -mtune=native")
  endif()
  set(CMAKE_CXX_FLAGS "${tmp}"
    CACHE STRING "cvt build defaults"
    FORCE)
  endif()
endif()

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

# VTK
# for VTK 8.2.0 use "vtk" and for 9.0.0 on use ""
set(vtk_ver 8)
set(vtk_mod_prefix)
if (vtk_ver EQUAL 8)
    set(vtk_mod_prefix vtk)
endif()
find_package(VTK COMPONENTS
    ${vtk_mod_prefix}CommonCore
    ${vtk_mod_prefix}CommonDataModel
    ${vtk_mod_prefix}IOCore
    ${vtk_mod_prefix}IOLegacy
    ${vtk_mod_prefix}IOParallelXML
    REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

# HDF5
find_package(HDF5 COMPONENTS C HL REQUIRED)
include_directories(${HDF5_C_INCLUDE_DIRS})

# Intel TBB
find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)

add_executable(cvt cvt.cpp neuron.cpp neuron_types.cpp)
target_link_libraries(cvt ${TBB_LIBRARIES}
    ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${VTK_LIBRARIES})

if (vtk_ver GREATER 8)
    vtk_module_autoinit(TARGETS cvt MODULES ${VTK_LIBRARIES})
endif()
