cmake_minimum_required(VERSION 3.4...3.18)
include(CMakePrintHelpers)
project(pyket LANGUAGES CXX)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_BUILD_TYPE Debug)
set(DEFAULT_NUM_MODES 2)
set(DEFAULT_NUM_BITS  8)

set(rapidjson_DIR ${CMAKE_CURRENT_SOURCE_DIRECTORY}/extern/rapidjson)
set(Eigen3_DIR /home/ethan/libraries/builds/eigen_3_install)
FIND_PACKAGE( Boost  )
FIND_PACKAGE( RapidJSON  )
find_package (Eigen3 3.3 NO_MODULE)



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

#find_path (${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11)
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt)
       execute_process(COMMAND git submodule update --init -- external/pybindll
                           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

endif()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/pyket_config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/pyket_config.hpp @ONLY)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11)	
pybind11_add_module(${pyket} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)	


if(NOT Boost_FOUND)
       message("unable to find boost, downloading instead")
       execute_process(COMMAND git submodule update --init -- external/boost-cmake
                           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	add_subdirectoy(${CMAKE_CURRENT_SOURCE_DIR}/external/boost-cmake)
       target_link_libraries( HamiltonianFuncs PRIVATE Boost::boost)
endif()

if(NOT rapidjson_FOUND )
       message("unable to find rapidjson, downloading instead")
       execute_process(COMMAND git submodule update --init extern/rapidjson
                           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	set(rapidJSON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/rapidjson CACHE PATH
	"rapidjson directory")
endif()

if(NOT Eigen3_FOUND)
       message("unable to find eigen3, downloading instead")
       execute_process(COMMAND git submodule update --init extern/eigen
                           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	set(Eigen3_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen/cmake CACHE PATH
	"Eigen3 directory")
       find_package (Eigen3 3.3 REQUIRED NO_MODULE)
endif()
target_link_libraries (HamiltonianFuncs Eigen3::Eigen)
target_compile_options(HamiltonianFuncs PRIVATE -fPIC)
target_link_libraries(${pyket} PRIVATE HamiltonianFuncs)
install(TARGETS ${pyket} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pyket)

