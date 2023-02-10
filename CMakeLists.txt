cmake_minimum_required(VERSION 3.5)
project(SnD VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(external/Argz)
add_subdirectory(external/wRATter)
add_subdirectory(external/spdlog)
add_subdirectory(external/nlopt)
add_subdirectory(external/boost)

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
        ${PROJECT_SOURCE_DIR}/external/wRATter/include
        ${PROJECT_SOURCE_DIR}/external/Argz/include
        ${PROJECT_SOURCE_DIR}/external/ProgressBar/include
        ${PROJECT_SOURCE_DIR}/external/ROOT/include
        ${ROOT_INCLUDE_DIRS}
)


#----------------------------------------------------------------------------
# LIBRARIES
#
add_library(Hit src/SnD/Hit.cc)
target_link_libraries(Hit ROOT::Physics)

add_library(NTuple src/Readers/NTuple.cc)
target_link_libraries(NTuple Hit ROOT::RIO ROOT::Tree Boost::range Boost::any)

add_library(Rat src/Readers/Rat.cc)
target_link_libraries(Rat Hit wRATter)


#----------------------------------------------------------------------------
# EXECUTABLES
#
add_executable(TApp src/Apps/TApp.cc src/Analyzers/TAnalyzer.cc)
target_link_libraries(TApp Argz NTuple ${ROOT_LIBRARIES})

add_executable(CreatePDF src/Apps/CreatePDF.cc src/Analyzers/MakePDF.cc)
target_link_libraries(CreatePDF Argz NTuple ${ROOT_LIBRARIES})

add_executable(Recon src/Apps/Recon.cc src/Analyzers/Recon.cc
        src/SnD/PosT.cc src/SnD/Multilateration.cc src/SnD/NLL.cc src/SnD/Recon.cc src/SnD/Map.cc src/SnD/SVD.cc)
target_link_libraries(Recon Argz NTuple ${ROOT_LIBRARIES} Boost::foreach Boost::optional nlopt)

