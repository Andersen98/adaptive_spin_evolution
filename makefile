objects = main.o input_tools.o
flags = -std=c++17 -g 
#-Wl,-t
adaptive_spin: ${objects}
	g++ ${flags}  -o adaptive_spin ${objects} -lboost_program_options  
main.o: configuration.hpp input_tools.hpp main.cpp hamiltonian.hpp spin_density_matrix.hpp
	g++ ${flags} -c main.cpp
input_tools.o: input_tools.hpp input_tools.cpp
	g++ ${flags} -c input_tools.cpp 

clean:
	rm -f *.o adaptive_spin
