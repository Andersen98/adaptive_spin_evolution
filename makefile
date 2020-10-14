objects = main.o input_tools.o
flags = -std=c++17

adaptive_spin: ${objects}
	g++ ${flags} -g -o adaptive_spin ${objects} -lboost_program_options
main.o: configuration.hpp input_tools.hpp main.cpp hamiltonian.hpp
	g++ ${flags} -c main.cpp
input_tools.o: input_tools.hpp input_tools.cpp
	g++ ${flags} -c input_tools.cpp 

clean:
	rm -f *.o adaptive_spin
