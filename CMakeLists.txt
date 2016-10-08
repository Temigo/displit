cmake_minimum_required(VERSION 2.8)
project(p3a)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

add_executable(main main.cpp)
target_link_libraries(main ${ROOT_LIBRARIES})