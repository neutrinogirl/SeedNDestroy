cmake_minimum_required(VERSION 3.5)
project(SeedNDestroy)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

## You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
##   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
###   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)
if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
endif()

find_package(Boost COMPONENTS system filesystem REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# ## NOT NEEDED FOR DICT
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/wRATter/include
        ${ROOT_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        $ENV{RATROOT}/include
        /usr/include/c++/4.9
        /usr/include/x86_64-linux-gnu/c++/4.9
)

#----------------------------------------------------------------------------
# Locate RATROOT Libs
#
LINK_DIRECTORIES($ENV{RATROOT}/lib)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(CreatePDF CreatePDF.cc)
add_executable(DebugRecon DebugRecon.cc)
add_executable(Recon Recon.cc)
