add_library(HamiltonianFuncs	STATIC
				  append_connections.cpp
				  merge_states.cpp
				  output.cpp
				  runtime.cpp
				  core.cpp
				  grow_configuration_space.cpp				   setup.cpp
				  input_tools.cpp
				  output_tools.cpp
				  )

target_compile_options(HamiltonianFuncs PRIVATE -fPIC)
target_include_directories(HamiltonianFuncs PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
