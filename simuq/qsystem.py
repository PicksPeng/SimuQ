from simuq.environment import BaseQuantumEnvironment

class QSystem(BaseQuantumEnvironment) :
    def __init__(self) :
        super().__init__()
        self.evos = []

    def add_evolution(self, h, t) :
        self.evos.append((h, t))



