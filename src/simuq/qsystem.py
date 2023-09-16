"""
The class for target Hamiltonian.

It is a container of the desired piecewise
constant evolution.
"""

from simuq.environment import BaseQuantumEnvironment, Boson, Fermion, Qubit


class QSystem(BaseQuantumEnvironment):
    """A target quantum system.

    It contains a list of evolutions (h, t), which
    represents evolving under h for time duration t.

    We also provide a syntax sugar to discretize
    a continuous-time Hamiltonian.
    """

    def __init__(self):
        super().__init__()
        self.evos = []

    def add_evolution(self, h, t):
        h.extend_ham_by_sites()
        self.evos.append((h, t))

    def add_time_dependent_evolution(self, ht, ts):
        for i in range(len(ts) - 1):
            self.add_evolution(ht(ts[i]), ts[i + 1] - ts[i])

    def add_td_evolution(self, ht, ts):
        return self.add_time_dependent_evolution(ht, ts)

    def clear_evos(self):
        self.evos = []

    def to_qutip(self):
        def time_indicator(tl, tr):
            def ret(t, args):
                if tl <= t and t < tr:
                    return 1
                else:
                    return 0

            return ret

        ret = []
        sumt = 0
        for h, t in self.evos:
            h.extend_ham_by_sites()
            ret.append([h.to_qutip_qobj(), time_indicator(sumt, sumt + t)])
            sumt += t

        return ret

    def total_time(self):
        ret = 0
        for h, t in self.evos:
            ret += t
        return ret

    def print_sites(self):
        name_list = []
        for site in self.sites:
            name_list.append(site.name)
        return name_list


if __name__ == "__main__":
    # This is to keep the classes under qsystem.
    tmp = Qubit(), Boson(), Fermion()
