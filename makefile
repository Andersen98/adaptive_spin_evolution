objects = main.o input_tools.o
flags = -std=c++17 -g 
dflags = -DNUM_MODES=1000 -DNUM_BITS=4

#-Wl,-t
adaptive_spin: ${objects}
	g++ ${dflags} ${flags}  -o adaptive_spin ${objects} -lboost_program_options  
main.o: configuration.hpp input_tools.hpp main.cpp hamiltonian.hpp spin_density_matrix.hpp
	g++ ${dflags} ${flags} -c main.cpp
input_tools.o: input_tools.hpp input_tools.cpp
	g++ ${dflags} ${flags} -c input_tools.cpp 

clean:
	rm -f *.o adaptive_spin
