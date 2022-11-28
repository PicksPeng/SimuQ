from simuq.qsystem import QSystem
from simuq.environment import qubit, fock

# Jaynes-Cummings model
qs = QSystem()
atom = qubit(qs)
cav = fock(qs)

omega_a = 0.1
omega_c = 0.01
g = 0.001
h = omega_a * atom.Z / 2 + omega_c * cav.c * cav.a + g * (cav.a + cav.c) * atom.X
t = 0.1
qs.add_evolution(h, t)
