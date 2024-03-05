from simuq.qmachine import QMachine
from simuq.environment import Qubit


def generate_qmachine(n=3):
    Ising = QMachine()
    qs = [Qubit(Ising) for _ in range(n)]
    h = [Ising.add_global_variable() for _ in range(n)]
    J = {(j, k): Ising.add_global_variable() for j in range(n) for k in range(n) if j > k}
    L = Ising.add_signal_line()
    ins = L.add_instruction("native")
    s = ins.add_local_variable()
    kinetic = s * sum([q.X for q in qs])
    problem = s * (sum([h[j] * qs[j].Z for j in range(n)]) + sum(
        [J[(j, k)] * qs[j].Z * qs[k].Z for j in range(n) for k in range(n) if
         j > k]))

    ins.set_ham(kinetic + problem)

    return Ising


generate_qmachine()
