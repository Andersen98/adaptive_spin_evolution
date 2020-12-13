objects = main.o input_tools.o output_tools.o
flags = -std=c++17 -g -O0
dflags = -DNUM_MODES=100 -DNUM_BITS=4

#-Wl,-t
adaptive_spin: ${objects}
	g++ ${dflags} ${flags}  -o adaptive_spin ${objects} -lboost_program_options  
main.o: configuration.hpp io_tools/input_tools.hpp io_tools/output_tools.hpp main.cpp hamiltonian.hpp spin_density_matrix.hpp
	g++ ${dflags} ${flags} -c main.cpp
input_tools.o: io_tools/input_tools.hpp io_tools/input_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/input_tools.cpp
output_tools.o: io_tools/output_tools.hpp io_tools/output_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/output_tools.cpp


clean:
	rm -f *.o adaptive_spin
