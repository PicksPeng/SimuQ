from simuq.environment import Qubit
from simuq.expression import Expression
from simuq.qmachine import QMachine


def generate_qmachine(n=3):
    mach = QMachine()
    ql = [Qubit(mach) for i in range(n)]
    link = [(i, i + 1) for i in range(n - 1)]

    for i in range(n):
        L = mach.add_signal_line()

        ins1 = L.add_instruction("native", "L{}_X_Y".format(i))
        amp = ins1.add_local_variable()
        phase = ins1.add_local_variable()
        ins1.set_ham(amp * (Expression.cos(phase) * ql[i].X + Expression.sin(phase) * ql[i].Y))

        ins2 = L.add_instruction("derived", "L{}_Z".format(i))
        amp = ins2.add_local_variable()
        ins2.set_ham(amp * ql[i].Z)

    for q0, q1 in link:
        L = mach.add_signal_line()

        ins = L.add_instruction("derived", "L{}{}_YX".format(q0, q1))
        amp = ins.add_local_variable()
        ins.set_ham(amp * ql[q0].Y * ql[q1].X)

    L = mach.add_signal_line()
    ins = L.add_instruction("derived", "L{}{}_YZX".format(q0, q1))
    amp = ins.add_local_variable()
    ins.set_ham(amp * ql[0].Y * ql[1].Z * ql[2].X)

    return mach
