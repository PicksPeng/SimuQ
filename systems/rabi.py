from simuq.qsystem import QSystem
from simuq.environment import qubit, fock
from simuq.qmachine import *
import math

qs = QSystem()
q = qubit(qs)
n_step = 5
T = 1.
omega0 = 0.1
omega1 = 0.3
omega = 0.2
for step in range(n_step) :
    t = step * T / n_step
    h = -omega0 / 2 * q.Z - omega1 / 2 * (math.cos(omega * t) * q.X - math.sin(omega * t) * q.Y)
    qs.add_evolution(h, T / n_step)
