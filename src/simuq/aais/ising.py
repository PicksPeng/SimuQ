from simuq.qmachine import QMachine
from simuq.environment import Qubit


def generate_qmachine(n=3):
    Ising = QMachine()
    qs = [Qubit(Ising) for _ in range(n)]
    L1 = Ising.add_signal_line()
    ins = L1.add_instruction("native")
    A = ins.add_local_variable()
    kinetic = (A / 2) * sum([q.X for q in qs])
    ins.set_ham(kinetic)
    h = [Ising.add_global_variable() for _ in range(n)]
    L2 = Ising.add_signal_line()
    ins = L2.add_instruction("native")
    B = ins.add_local_variable()
    problem = (B / 2) * (sum([h[j] * qs[j].Z for j in range(n)]) + sum(
        [Ising.add_global_variable() * qs[j].Z * qs[k].Z for j in range(n) for k in range(n) if
         j > k]))
    ins.set_ham(problem)
    return Ising


generate_qmachine()
