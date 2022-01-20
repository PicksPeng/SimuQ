from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression

Rydberg = QMachine()

q1 = qubit(Rydberg)
q2 = qubit(Rydberg)
q3 = qubit(Rydberg)

x0 = 0
x1 = GlobalVar(Rydberg)
x2 = GlobalVar(Rydberg)

L1 = SignalLine(Rydberg)

ins1 = Instruction(L1, 'native', 'ins1')
a = LocalVar(ins1)
ins1.set_ham((Expression.cos(a) + 2 / (x2 - x1) ** 4) * q1.X() * q2.X())

ins2 = Instruction(L1, 'derived', 'ins2')
a = LocalVar(ins2)
ins2.set_ham((a * (x1 + x0)) * (q2.X() + q1.Y()))
