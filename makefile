objects = main.o input_tools.o output_tools.o
flags = -fPIC -std=c++20 -ggdb3 -O0 -Wall 
dflags = -DNUM_MODES=1000 -DNUM_BITS=8
hamiltonian_objects = grow_configuration_space.o evolve_space.o append_connections.o merge_states.o setup.o core.o output.o
objects += ${hamiltonian_objects}
#-Wl,-t
pyket_objects += ${hamiltonian_objects} pyket.o input_tools.o
pyket_flags += ${dflags} ${flags}  -shared -I extern/pybind11/include `python-config --includes`
pyket_name = pyket`python3-config --extension-suffix`

pyket: ${pyket_objects}
	g++ ${pyket_flags} -o job_manager/${pyket_name} ${pyket_objects} -lboost_program_options
pyket.o: pyket/pyket.cpp hamiltonian/hamiltonian.hpp configuration.hpp
	g++ ${pyket_flags} -c pyket/pyket.cpp


adaptive_spin: ${objects}
	g++ ${dflags} ${flags}  -o adaptive_spin ${objects} -lboost_program_options  
main.o: configuration.hpp io_tools/input_tools.hpp io_tools/output_tools.hpp main.cpp hamiltonian/hamiltonian.hpp spin_density_matrix.hpp
	g++ ${dflags} ${flags} -c main.cpp
input_tools.o: io_tools/input_tools.hpp io_tools/input_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/input_tools.cpp
output_tools.o:configuration.hpp io_tools/output_tools.hpp io_tools/output_tools.cpp
	g++ ${dflags} ${flags} -c io_tools/output_tools.cpp

grow_configuration_space.o: configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/grow_configuration_space.cpp
	g++ ${dflags} ${flags} -c hamiltonian/grow_configuration_space.cpp

append_connections.o:configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/append_connections.cpp
	g++ ${dflags} ${flags} -c hamiltonian/append_connections.cpp

evolve_space.o: configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/evolve_space.cpp
	g++ ${dflags} ${flags} -c hamiltonian/evolve_space.cpp

merge_states.o: configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/merge_states.cpp
	g++ ${dflags} ${flags} -c hamiltonian/merge_states.cpp

setup.o:configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/setup.cpp io_tools/input_tools.hpp
	g++ ${dflags} ${flags} -c hamiltonian/setup.cpp

core.o: configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/core.cpp
	g++ ${dflags} ${flags} -c hamiltonian/core.cpp

output.o:configuration.hpp hamiltonian/hamiltonian.hpp hamiltonian/output.cpp
	g++ ${dflags} ${flags} -c hamiltonian/output.cpp

clean:
	rm -f *.o adaptive_spin job_manager/${pyket_name}


