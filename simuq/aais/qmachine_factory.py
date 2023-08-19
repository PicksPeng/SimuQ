from abc import ABC, abstractmethod
from simuq.qmachine import QMachine

class QMachineFactory(ABC):
    @staticmethod
    @abstractmethod
    def generate_qmachine(*args, **kwargs) -> QMachine:
        pass

