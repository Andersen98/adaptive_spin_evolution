cmake_minimum_required(VERSION 3.4...3.18)
project(adaptiv_evolution LANGUAGES CXX)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(rapidjson_DIR /home/ethan/libraries/builds/rapidjson)

FIND_PACKAGE( Boost  REQUIRED )
FIND_PACKAGE( RapidJSON  REQUIRED )

add_library(HamiltonianFuncs "")
add_subdirectory(src)


add_subdirectory(pybind11)

pybind11_add_module(pyket src/main.cpp)

target_link_libraries(pyket PUBLIC HamiltonianFuncs)
target_compile_definitions(pyket PUBLIC NUM_MODES=1000 NUM_BITS=8)

