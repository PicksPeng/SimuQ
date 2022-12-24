'''
The class for target Hamiltonian.

It is a container of the desired piecewise
constant evolution.
'''

from simuq.environment import BaseQuantumEnvironment

class QSystem(BaseQuantumEnvironment) :
    """ A target quantum system.

    It contains a list of evolutions (h, t), which 
    represents evolving under h for time duration t.

    We also provide a syntax sugar to discretize
    a continuous-time Hamiltonian.
    """
    def __init__(self) :
        super().__init__()
        self.evos = []

    def add_evolution(self, h, t) :
        self.evos.append((h, t))

    def add_time_dependent_evolution(self, ht, ts) :
        for i in range(len(ts) - 1) :
            self.add_evolution(ht(ts[i]), ts[i + 1] - ts[i])

    def clear_evos(self) :
        self.evos = []
