add_library(IO_Tools
		input_tools.cpp
		output_tools.cpp
		)
target_include_directories(IO_Tools
		PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	       	)

# install rules
install(TARGETS IO_Tools DESTINATION lib)
install(FILES input_tools.hpp output_tools.hpp DESTINATION include)
