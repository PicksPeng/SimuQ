from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression

mach = QMachine()
n = 7
ql = [qubit(mach) for i in range(27)]
link = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 9), (8, 11), (11, 14), (1, 4), (4, 7), (6, 7), (7, 10), (10, 12), (12, 13), (13, 14), (12, 15), (15, 18), (17, 18), (18, 21), (21, 23), (23, 24), (14, 16), (16, 19), (19, 20), (19, 22), (22, 25), (24, 25), (25, 26)]

for i in range(n) :
    L = SignalLine(mach)
    
    ins1 = Instruction(L, 'native', 'L{}_X_Y'.format(i))
    amp = LocalVar(ins1)
    phase = LocalVar(ins1)
    ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))
    
    ins2 = Instruction(L, 'derived', 'L{}_Z'.format(i))
    amp = LocalVar(ins2)
    ins2.set_ham(amp * ql[i].Z)

for (q0, q1) in link :
    L = SignalLine(mach)

    ins = Instruction(L, 'derived', 'L{}{}_ZX'.format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Z * ql[q1].X)

    ins = Instruction(L, 'derived', 'L{}{}_XX'.format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].X * ql[q1].X)

    ins = Instruction(L, 'derived', 'L{}{}_YY'.format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Y * ql[q1].Y)

    ins = Instruction(L, 'derived', 'L{}{}_ZZ'.format(q0, q1))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].Z * ql[q1].Z)

    L = SignalLine(mach)
    
    ins = Instruction(L, 'derived', 'L{}{}_ZX'.format(q1, q0))
    amp = LocalVar(ins)
    ins.set_ham(amp * ql[q0].X * ql[q1].Z)

