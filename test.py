# from aais.rydberg1d import Rydberg as mach
from aais.ibm import mach
from backends.bloqade_rydberg import transpile
from simuq.solver import generate_as
from systems.qaoa import qs

# from systems.mis import qs
# from aais.iontrap import mach
# from backends.qiskit_iontrap import transpile
# from aais.ibm5 import mach
# from backends.qiskit_5 import transpile


# from systems.heis5 import qs
# from systems.single_x import qs
# from systems.annealing import qs


# generate_as(FloquetQS, FluxSCMach)
# generate_as(qs, mach, 1)


# from aais.rydberg2d import Rydberg as mach
# from backends.bloqade_rydberg2d import transpile


# generate_as(qs, mach, 1, tol = 1e-1)

print(generate_as(qs, mach, 1))

"""
circ=transpile(*generate_as(qs, mach, 1, 0.5))
import pickle

with open('circ.pickle', 'wb') as f:
    pickle.dump(circ, f)
"""
