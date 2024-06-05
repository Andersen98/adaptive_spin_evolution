objects = input_tools.o main_test.o
flags = -g -std=c++17
#-DBOOST_TEST_TOOLS_UNDER_DEBUGGER
boost_test: ${objects}
	g++ ${flags} -o boost_test ${objects} -lboost_unit_test_framework -lboost_program_options

main_test.o: ../configuration.hpp ../input_tools.hpp ../hamiltonian.hpp main_test.cpp
	g++ ${flags} -c main_test.cpp

input_tools.o: ../input_tools.hpp ../input_tools.cpp
	g++ ${flags} -c ../input_tools.cpp
clean:
	rm *.o boost_test


