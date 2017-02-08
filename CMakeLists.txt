cmake_minimum_required(VERSION 2.8)
project(p3a)

# Use MPI
SET(CMAKE_CXX_COMPILER mpicxx)
# Use C++11
set (CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT COMPONENTS MathMore REQUIRED)
include(${ROOT_USE_FILE})

#find_package(GSL REQUIRED)
#include_directories(${GSL_INCLUDE_DIRS})

add_library(globals globals.cpp)
target_include_directories(globals PUBLIC include/ )

add_library(EventTree EventTree.C)
target_link_libraries(EventTree ${ROOT_LIBRARIES})
target_include_directories(EventTree PUBLIC include/ )

add_library(dipole dipole.cpp)
target_link_libraries(dipole ${ROOT_LIBRARIES})
target_include_directories(dipole PUBLIC include/ )

add_library(event event.cpp)
target_link_libraries(event ${ROOT_LIBRARIES} dipole globals)
target_include_directories(event PUBLIC include/ )

add_library(utils utils.cpp)
target_link_libraries(utils ${ROOT_LIBRARIES} event EventTree globals)
target_include_directories(utils PUBLIC include/ )

add_library(graphics graphics.cpp)
target_link_libraries(graphics ${ROOT_LIBRARIES} event)
target_include_directories(graphics PUBLIC include/ )

add_executable(main main.cpp)
target_link_libraries(main ${ROOT_LIBRARIES} utils graphics mpi)
target_include_directories(main PUBLIC include/ )