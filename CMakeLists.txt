cmake_minimum_required(VERSION 3.11)
project(SnD VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(external/Argz)
add_subdirectory(external/range-v3)
add_subdirectory(external/nlopt)
add_subdirectory(external/googletest)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT 6.24 CONFIG REQUIRED)
if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# ## NOT NEEDED FOR DICT
include_directories(
        ${PROJECT_SOURCE_DIR}/include
#        ${PROJECT_SOURCE_DIR}/external/wRATter/include
        ${PROJECT_SOURCE_DIR}/external/Argz/include
        ${PROJECT_SOURCE_DIR}/external/ProgressBar/include
        ${PROJECT_SOURCE_DIR}/external/ROOT/include
        ${PROJECT_SOURCE_DIR}/external/range-v3/include
        ${ROOT_INCLUDE_DIRS}
)


#----------------------------------------------------------------------------
# LIBRARIES
#

add_library(SVD src/LinAlg/SVD.cc)

add_library(Hit src/SnD/Hit.cc)
target_link_libraries(Hit SVD ROOT::Hist)

add_library(TReader src/Templates/TReader.cc)

add_library(NTuple src/Readers/NTuple.cc)
target_link_libraries(NTuple TReader Hit ROOT::RIO ROOT::Tree ROOT::TreePlayer)
#
#add_library(Rat src/Readers/Rat.cc)
#target_link_libraries(Rat Hit wRATter)


#----------------------------------------------------------------------------
# EXECUTABLES
#
add_executable(TApp src/Apps/TApp.cc src/Analyzers/TAnalyzer.cc)
target_link_libraries(TApp Argz NTuple)

add_executable(CreatePDF src/Apps/CreatePDF.cc src/Analyzers/MakePDF.cc)
target_link_libraries(CreatePDF Argz NTuple)

add_executable(Recon src/Apps/Recon.cc src/Analyzers/Recon.cc src/Algo/VHits.cc src/Algo/WOpt.cc)
target_link_libraries(Recon Argz NTuple nlopt)

#add_executable(Map src/Apps/Map.cc src/Analyzers/Map.cc
#        src/Algo/VHits.cc
#        src/Algo/WOpt.cc
#        src/SnD/SpaceGrid.cc)
#target_link_libraries(Map Argz NTuple nlopt)

#----------------------------------------------------------------------------
# TESTS
#
add_executable(Vector3Test test/Vector3Test.cc)
target_link_libraries(Vector3Test gtest_main)

add_executable(HitTest test/HitTest.cc)
target_link_libraries(HitTest gtest_main Hit)