# we need at least cmake version 3.1
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(qmc_worm CXX)

add_compile_options(
  "-Wall" "-Wextra"
  "-Wno-unused-parameter"
  "-Wno-reorder"
  "-Wfatal-errors"
  "-Wno-unknown-pragmas"
  "-fdiagnostics-color=always"
  )

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(DEBUG OFF CACHE BOOL "Build in Debug Mode")
if(DEBUG)
  add_definitions(-DDEBUGMODE)
  set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "Build type is Debug")
else()
  remove_definitions(-DDEBUGMODE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type is Release")
endif()

set(UNIFORM ON CACHE BOOL "System is uniform")
if (UNIFORM)
  add_definitions(-DUNISYS)
  message(STATUS "System is uniform")
else()
  remove_definitions(-DUNISYS)
  message(STATUS "System is not uniform")
endif()

set(DSTRUC LIST CACHE STRING "Data Structure")
if(DSTRUC STREQUAL AVL)
  add_definitions(-DDSTRUC_AVL)
  remove_definitions(-DDSTRUC_LIST)
  remove_definitions(-DDSTRUC_LIST_STACK)
  message(STATUS "Using Data Structure AVL-tree")
elseif(DSTRUC STREQUAL LIST)
  remove_definitions(-DDSTRUC_AVL)
  add_definitions(-DDSTRUC_LIST)
  remove_definitions(-DDSTRUC_LIST_STACK)
  message(STATUS "Using Data Structure LIST")
elseif(DSTRUC STREQUAL LIST_STACK)
  remove_definitions(-DDSTRUC_AVL)
  remove_definitions(-DDSTRUC_LIST)
  add_definitions(-DDSTRUC_LIST_STACK)
  message(STATUS "Using Data Structure LIST_STACK")
endif()

set(CWINDOW OFF CACHE BOOL "Canonical Window for Green's Function")
if(CWINDOW)
    add_definitions(-DCAN_WINDOW)
    message(STATUS "Canonical Window for GF computation is set")
else()
    remove_definitions(-DCAN_WINDOW)
    message(STATUS "Canonical Window for GF computation is not set")
endif()


add_executable(${PROJECT_NAME}
  worm.run.cpp 
  worm.cpp 
  worm.update.cpp 
  worm.output.cpp)

add_executable(${PROJECT_NAME}_mpi
  worm.run_mpi.cpp
  worm.cpp
  worm.update.cpp
  worm.output.cpp)


# the project relies on the ALPSCore package. If ALPSCore is not found
# automatically, specify its location using:
# export ALPSCore_DIR=/location/to/ALPSCORE/
find_package(ALPSCore 2.0  REQUIRED mc accumulators params)

# check if Eigen is present:
if (ALPSCore_HAS_EIGEN_VERSION)
  message(STATUS "ALPSCore provides Eigen, version ${ALPSCore_HAS_EIGEN_VERSION}")
else()
  message("WARNING: ALPSCore does not provide Eigen, the code won't compile!")
endif()
# Use ALPSCore_LIBRARIES variable to link to ALPSCore 

target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})
target_link_libraries(${PROJECT_NAME}_mpi ${ALPSCore_LIBRARIES})
