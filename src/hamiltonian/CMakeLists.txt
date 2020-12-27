target_sources(HamiltonianFuncs  PRIVATE
				  append_connections.cpp
				  evolve_space.cpp
				  merge_states.cpp
				  output.cpp
				  runtime.cpp
				  core.cpp
				  grow_configuration_space.cpp
				  odeint_evolve.cpp
				  par_runtime.cpp
				  setup.cpp
				  io_tools/input_tools.hpp
				  io_tools/input_tools.cpp
				  PUBLIC
				  configuration.hpp
				  )
				  
target_compile_definitions(HamiltonianFuncs PUBLIC NUM_MODES=1000 NUM_BITS=8)
target_compile_options(HamiltonianFuncs PUBLIC -fPIC )
target_include_directories(HamiltonianFuncs PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/io_tools )

