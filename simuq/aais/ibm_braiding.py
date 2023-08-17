from simuq.aais.qmachine_factory import QMachineFactory
from simuq.qmachine import *

class IbmBraidingFactory(QMachineFactory):
    @staticmethod
    def generate_qmachine(n=3, *args, **kwargs):
        mach = QMachine()
        n = 3
        ql = [qubit(mach) for i in range(n)]
        link = [(i, i + 1) for i in range(n - 1)]

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

            ins = Instruction(L, 'derived', 'L{}{}_YX'.format(q0, q1))
            amp = LocalVar(ins)
            ins.set_ham(amp * ql[q0].Y * ql[q1].X)

        L = SignalLine(mach)
        ins = Instruction(L, 'derived', 'L{}{}_YZX'.format(q0, q1))
        amp = LocalVar(ins)
        ins.set_ham(amp * ql[0].Y * ql[1].Z * ql[2].X)

        return mach