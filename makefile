objects = main.o input_tools.o output_tools.o
flags = -std=c++20 -ggdb3 -O0 -Wall
dflags = -DNUM_MODES=2 -DNUM_BITS=8
hamiltonian_objects = grow_configuration_space.o evolve_space.o append_connections.o merge_states.o setup.o core.o output.o
objects += ${hamiltonian_objects}
#-Wl,-t
adaptive_spin: ${objects}
	g++ ${dflags} ${flags}  -o adaptive_spin ${objects} -lboost_program_options  
main.o: configuration.hpp io_tools/input_tools.hpp io_tools/output_tools.hpp main.cpp hamiltonian/hamiltonian.hpp spin_density_matrix.hpp
	g++ ${dflags} ${flags} -c main.cpp
input_tools.o: io_tools/input_tools.hpp io_tools/input_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/input_tools.cpp
output_tools.o: io_tools/output_tools.hpp io_tools/output_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/output_tools.cpp

grow_configuration_space.o: hamiltonian/hamiltonian.hpp hamiltonian/append_connections.cpp
	g++ ${dflags} ${flags} -c hamiltonian/grow_configuration_space.cpp

append_connections.o: hamiltonian/hamiltonian.hpp hamiltonian/append_connections.cpp
	g++ ${dflags} ${flags} -c hamiltonian/append_connections.cpp

evolve_space.o: hamiltonian/hamiltonian.hpp hamiltonian/evolve_space.cpp
	g++ ${dflags} ${flags} -c hamiltonian/evolve_space.cpp

merge_states.o: hamiltonian/hamiltonian.hpp hamiltonian/merge_states.cpp
	g++ ${dflags} ${flags} -c hamiltonian/merge_states.cpp

setup.o: hamiltonian/hamiltonian.hpp hamiltonian/setup.cpp
	g++ ${dflags} ${flags} -c hamiltonian/setup.cpp

core.o: hamiltonian/hamiltonian.hpp hamiltonian/core.cpp
	g++ ${dflags} ${flags} -c hamiltonian/core.cpp

output.o: hamiltonian/hamiltonian.hpp hamiltonian/output.cpp
	g++ ${dflags} ${flags} -c hamiltonian/output.cpp

clean:
	rm -f *.o adaptive_spin
