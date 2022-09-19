from simuq.solver import generate_as

#from systems.heis5 import qs
#from systems.single_x import qs
#from systems.annealing import qs

#from aais.ibm7 import mach

#generate_as(FloquetQS, FluxSCMach)
#generate_as(qs, mach, 1)



from systems.mis import qs
from aais.rydberg1d import Rydberg as mach
#from aais.rydberg2d import Rydberg as mach
from backends.bloqade_rydberg import transpile
#from backends.bloqade_rydberg2d import transpile


#from systems.mis import qs
#from aais.iontrap import mach
#from backends.qiskit_iontrap import transpile

#generate_as(qs, mach, 1, tol = 1e-1)

print(transpile(*generate_as(qs, mach, 1, 0.5)))
