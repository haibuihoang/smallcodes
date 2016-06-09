FC = g95
sevol.exe:  arat2p.o arat2v.o se.o comm.mod data.mod  module_solve.o plot.mod
	${FC} -o sevol.exe arat2p.o arat2v.o module_comm.o module_data.o module_solve.o se.o module_plot.o

arat2p.o: arat2p.f90
	${FC} -c arat2p.f90 
arat2v.o: arat2v.f90
	${FC} -c arat2v.f90 
comm.mod: module_comm.f90
	${FC} -c module_comm.f90 
data.mod: comm.mod module_data.f90
	${FC} -c module_data.f90 
plot.mod: data.mod solve.mod module_plot.f90
	${FC} -c module_plot.f90 
module_solve.o: data.mod comm.mod module_solve.f90
	${FC} -c module_solve.f90 
se.o: se.f90 data.mod module_solve.o plot.mod
	${FC} -c se.f90 

clean:
	rm *.exe *.o *.mod