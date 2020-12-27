cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
add_subdirectory(pybind11)
pybind11_add_module(pyket src/main.cpp)
project(adaptiv_evolution LANGUAGES CXX)

add_compile_definitions(NUM_MODES=1000 NUM_BITS=8)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(rapidjson_DIR /home/ethan/libraries/builds/rapidjson)
set(eigen_DIR /home/ethan/libraries/builds/eigen)
FIND_PACKAGE( Boost  REQUIRED )
FIND_PACKAGE( RapidJSON  REQUIRED )
FIND_PACKAGE( Eigen  REQUIRED )
add_library(HamiltonianFuncs SHARED "")
add_subdirectory(src)



target_link_libraries(pyket PUBLIC HamiltonianFuncs)


