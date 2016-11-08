cmake_minimum_required(VERSION 2.8)
project(p3a)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

add_library(MyClass MyClass.C)
target_link_libraries(MyClass ${ROOT_LIBRARIES})
target_include_directories(MyClass PUBLIC include/ )

add_library(dipole dipole.cpp)
target_link_libraries(dipole ${ROOT_LIBRARIES})
target_include_directories(dipole PUBLIC include/ )

add_library(event event.cpp)
target_link_libraries(event ${ROOT_LIBRARIES} dipole)
target_include_directories(event PUBLIC include/ )

add_executable(main main.cpp)
target_link_libraries(main ${ROOT_LIBRARIES} event MyClass)
target_include_directories(main PUBLIC include/ )