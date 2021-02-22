cmake_minimum_required(VERSION 3.6.1)
project(CiSE-Pandemic)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_COMPILER "g++")
add_compile_options(-g)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include_directories(cadmium/include)
include_directories(cadmium/json/include)

set(Boost_USE_MULTITHREADED TRUE)
find_package(Boost COMPONENTS unit_test_framework system thread REQUIRED)

file(MAKE_DIRECTORY simulation_results)

add_executable(CiSE-Pandemic model/main.cpp)

target_link_libraries(CiSE-Pandemic PUBLIC ${Boost_LIBRARIES})
