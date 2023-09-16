from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import *

mach = QMachine()
n = 7
ql = [Qubit(mach) for i in range(7)]
link = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]

for i in range(n):
    L = SignalLine(mach)

    ins1 = Instruction(L, "native", "L{}_X_Y".format(i))
    amp = LocalVar(ins1)
    phase = LocalVar(ins1)
    ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))

    ins2 = Instruction(L, "derived", "L{}_Z".format(i))
    amp = LocalVar(ins2)
    ins2.set_ham(amp * ql[i].Z)

for q0, q1 in link:
    L = SignalLine(mach)

    ins = Instruction(L, "derived", "L{}{}_ZX".format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Z * ql[q1].X)

    ins = Instruction(L, "derived", "L{}{}_XX".format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].X * ql[q1].X)

    ins = Instruction(L, "derived", "L{}{}_YY".format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Y * ql[q1].Y)

    ins = Instruction(L, "derived", "L{}{}_ZZ".format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Z * ql[q1].Z)

    L = SignalLine(mach)

    ins = Instruction(L, "derived", "L{}{}_ZX".format(q1, q0))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].X * ql[q1].Z)
