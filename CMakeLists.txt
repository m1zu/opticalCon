cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(dummyRootProject)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO)

#---Locate Qt
message(STATUS "Searching for Qt5")
FIND_PACKAGE(Qt5Core QUIET)
if (Qt5Core_FOUND)
  set(QT_FOUND true)
  add_definitions(${Qt5Core_DEFINITIONS})
  include_directories(${Qt5Core_INCLUDE_DIRS})
  set(QT_INCLUDES ${Qt5Core_INCLUDE_DIRS})
  set(QT_LIBRARIES ${Qt5Core_LIBRARIES})
  set(QT_VERSIONSTRING ${Qt5Core_VERSION_STRING})
  message(STATUS "Found Qt5 version ${QT_VERSIONSTRING}")
  set(QTV "5")
endif()
if (NOT QT_FOUND)
  message(STATUS "Qt5 not found")
  message(STATUS "Searching for Qt4")
  find_package(Qt4 COMPONENTS QtCore QUIET)
  if (QT4_FOUND)
    set(QT_FOUND true)
    add_definitions(${QT4_DEFINITIONS})
    include(${QT_USE_FILE})
    set(QT_VERSIONSTRING "${QT_VERSION_MAJOR}.${QT_VERSION_MINOR}.${QT_VERSION_PATCH}")
    set(QT_VERSIONSTRING "${QT_VERSIONSTRING}")
    message(STATUS "Found Qt4 version ${QT_VERSIONSTRING}")
    set(QTV "4")
  endif()
  if (NOT QT_FOUND)
    message(STATUS "Qt4 not found")
    message(FATAL_ERROR "Neither Qt4 nor Qt5 were found, but one of them is required.")
  endif()
endif()

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(source_files
  RootExecFunctionWrapper.cc
  Functions.cc
  Helpers.cc
  main.cc
)

set(header_files
  RootExecFunctionWrapper.hh
  Functions.hh
  Helpers.hh
)

#---Create a main program using the libraries
add_executable(Main ${source_files} ${header_files})
target_link_libraries(Main ${ROOT_LIBRARIES} ${QT_LIBRARIES})
