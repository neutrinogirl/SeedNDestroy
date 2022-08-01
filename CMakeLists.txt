cmake_minimum_required(VERSION 3.5)
project(SnD VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(external/wRATter)
add_subdirectory(external/Argz)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT 6.24 CONFIG REQUIRED)
if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
endif()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.50 REQUIRED COMPONENTS filesystem)
message(STATUS "Boost version: ${Boost_VERSION}")

set(NLOPT_INCLUDE_DIRS ${NLOPT_ROOT}/include)
set(NLOPT_LIB ${NLOPT_ROOT}/lib)

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
        ${Boost_INCLUDE_DIRS}
        ${NLOPT_INCLUDE_DIRS}
        $ENV{RATROOT}/include
)

#----------------------------------------------------------------------------
# Locate RATROOT Libs
#
LINK_DIRECTORIES($ENV{RATROOT}/lib)
LINK_DIRECTORIES($ENV{GEANT4SYS}/lib)
LINK_DIRECTORIES(${NLOPT_LIB})

###
# ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${NLOPT_LIB}

#----------------------------------------------------------------------------
# LIBRARIES
#
add_library(Hit src/SnD/Hit.cc)
target_link_libraries(Hit ROOT::Physics)

add_library(RATReader src/SnD/RATReader.cc)
target_link_libraries(RATReader wRATter)

add_library(SVD src/SnD/SVD.cc)

add_library(Multilateration src/SnD/Multilateration.cc)
target_link_libraries(Multilateration ROOT::Physics ROOT::Hist)

#----------------------------------------------------------------------------
# EXECUTABLES
#
add_executable(TApp src/CreatePDF/TApp.cc src/CreatePDF/PDFAnalysis.cc)
target_link_libraries(TApp Argz RATReader Hit ROOT::Physics ROOT::Hist ROOT::RIO ROOT::Tree)

add_executable(ReconApp src/Recon/ReconApp.cc src/Recon/ReconAnalysis.cc src/SnD/PosT.cc src/SnD/Multilateration.cc src/SnD/NLL.cc)
target_link_libraries(ReconApp Argz RATReader Hit SVD ROOT::Physics ROOT::Hist ROOT::RIO ROOT::Tree)

add_executable(ERec src/ERec/TApp.cc src/ERec/ERecAnalysis.cc)
target_link_libraries(ERec Argz RATReader Hit ROOT::Physics ROOT::Hist ROOT::RIO ROOT::Tree)
