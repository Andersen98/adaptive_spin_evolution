cmake_minimum_required(VERSION 3.4...3.18)
include(CMakePrintHelpers)
project(pyket LANGUAGES CXX)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
#set(CMAKE_BUILD_TYPE Debug)
set(DEFAULT_NUM_MODES 2)
set(DEFAULT_NUM_BITS  8)

add_subdirectory(extern/boost)
add_subdirectory(extern/eigen)

#set(rapidjson_DIR home/ethan/libraries/builds/rapidjson/include)
#set(Eigen3_DIR /home/ethan/libraries/builds/eigen_3_install)
#FIND_PACKAGE( Boost REQUIRED )
#find_package(Eigen3 3.3 NO_MODULE REQUIRED)



if(NOT DEFINED NUM_MODES)
       message(STATUS "Using Default Number of modes")
       set(NUM_MODES  ${DEFAULT_NUM_MODES})
endif()
if(NOT DEFINED NUM_BITS)
       message(STATUS "Using Default Number of bits per mode")
       set(NUM_BITS  ${DEFAULT_NUM_BITS})
endif()

#-------DEFINE MODULE NAME------------
set(pyket pyket${NUM_MODES}_${NUM_BITS})

add_compile_definitions(NUM_MODES=${NUM_MODES} NUM_BITS=${NUM_BITS})
cmake_print_variables(NUM_MODES)
cmake_print_variables(NUM_BITS)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_library(HamiltonianFuncs STATIC "")

add_subdirectory(src)


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11)	
pybind11_add_module(${pyket} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)	


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/pyket_config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/pyket_config.hpp @ONLY)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(HamiltonianFuncs SYSTEM PUBLIC extern/rapidjson/include)
target_link_libraries(HamiltonianFuncs Eigen3::Eigen)
target_link_libraries(HamiltonianFuncs Boost::format)
target_compile_options(HamiltonianFuncs PRIVATE -fPIC)
target_link_libraries(${pyket} PRIVATE HamiltonianFuncs)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${pyket} PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS ${pyket} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pyket)

