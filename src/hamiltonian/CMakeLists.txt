target_sources(HamiltonianFuncs  PRIVATE
				  append_connections.cpp
				  merge_states.cpp
				  output.cpp
				  runtime.cpp
				  core.cpp
				  grow_configuration_space.cpp				  
				  par_runtime.cpp
				  setup.cpp
				  io_tools/input_tools.hpp
				  io_tools/input_tools.cpp
				  PUBLIC
				  configuration.hpp
				  )
				  

target_include_directories(HamiltonianFuncs PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/io_tools )

