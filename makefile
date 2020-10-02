objects = main.o input_tools.o

adaptive_spin: ${objects}
	g++ -g -o adaptive_spin ${objects} -lboost_program_options
main.o: configuration.hpp input_tools.hpp main.cpp hamiltonian.hpp
	g++ -c main.cpp
input_tools.o: input_tools.hpp input_tools.cpp
	g++ -c input_tools.cpp 

clean:
	rm -f *.o adaptive_spin
