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
  set(tmp "-fPIC -std=c++11 -Wall -Wextra")
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

find_package(VTK COMPONENTS CommonCore CommonDataModel IOCore IOLegacy IOParallelXML REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)

add_executable(cvt cvt.cpp)
target_link_libraries(cvt ${TBB_LIBRARIES}
    ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${VTK_LIBRARIES})

vtk_module_autoinit(TARGETS cvt MODULES ${VTK_LIBRARIES})
