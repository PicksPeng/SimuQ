"""
The class file for quantum machines.

One can implement an AAIS in a QMachine. Signal lines
and instructions are supported inside it.

A quantum machine is essentially a container of signal
lines, instructions and (local or global) variables.

A signal line contains multiple instructions. Instructions
sharing a signal line cannot be called simultaneously.

Global variables are defined for the QMachine as a
globally tunable parameter. They will be fixed for the
whole experiment.

An instruction may contain several local variables. Local
variables can be tuned for each call of the instruction.
Both variables can be used in the description of its
Hamiltonian.
TODO: Add API for variables' bounds.
"""

from copy import copy

import numpy as np

from simuq.environment import BaseQuantumEnvironment
from simuq.expression import BaseVar, Expression


class QMachine(BaseQuantumEnvironment):
    """A quantum device described by its AAIS.

    It stores the signal lines, instructions, local
    variables and global variables.
    """

    def __init__(self):
        super().__init__()
        self.gvars = []
        self.lvars = []
        self.lines = []
        self.num_inss = 0
        self.sys_ham = 0
        self.with_sys_ham = False
        self.instantiated = False

    def set_sys_ham(self, h):
        self.sys_ham = h
        self.with_sys_ham = True

    def instantiate_sys_ham(self):
        if self.with_sys_ham and not self.instantiated:
            SysLine = SignalLine(self)
            SysIns = Instruction(SysLine, "native", "System Hamiltonian")
            SysIns.set_ham(self.sys_ham)
            SysIns.is_sys_ham = True
            self.instantiated = True

    def extend_instruction_sites(self):
        for line in self.lines:
            for ins in line.inss:
                ins.h.extend_sites(self.sites_type)

    def add_signal_line(self):
        return SignalLine(self)

    def add_global_variable(self, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        return GlobalVar(self, init_value, lower_bound, upper_bound)


class SignalLine:
    """A signal line.

    It contains all instructions belonging to it.
    """

    def __init__(self, mach):
        self.mach = mach
        self.index = len(mach.lines)
        mach.lines.append(self)
        self.inss = []

    def add_instruction(self, nativeness="native", name="ins"):
        return Instruction(self, nativeness, name)


class Instruction:
    """An instruction.

    It contains the local variables belonging to it, its
    property, and its instruction Hamiltonian.
    """

    def __init__(self, line, nativeness="native", name="ins"):
        self.line = line
        self.line.inss.append(self)
        self.index = self.line.mach.num_inss
        self.line.mach.num_inss += 1
        self.vars_index = []
        self.nativeness = nativeness
        self.name = name
        self.is_sys_ham = False
        self.h = None

    @property
    def mach(self):
        return self.line.mach

    def set_ham(self, h):
        newh = copy(h)
        for i in range(len(h.ham)):
            if not isinstance(h.ham[i][1], Expression):
                newh.ham[i] = (h.ham[i][0], h.ham[i][1] * Expression.unit(self.mach))
        self.h = newh

    def exp_eval(self, gvars, lvars):
        return self.h.exp_eval(gvars, lvars)

    def add_local_variable(self, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        return LocalVar(self, init_value, lower_bound, upper_bound)


class GlobalVar(BaseVar):
    """A global variable, belonging to a QMachine.

    One can set its initial value, lower and upper bounds
    when declaring it.
    """

    def __init__(self, mach, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        super().__init__(mach)
        self.index = len(mach.gvars)
        mach.gvars.append(self)
        self.init_value = init_value
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def to_exp(self):
        e = Expression.id_gvar(self.mach, self.index)
        return e


class LocalVar(BaseVar):
    """A local vabiable, belonging to an instruction.

    One can set its initial value, lower and upper bounds
    when declaring it.
    """

    def __init__(self, ins, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        mach = ins.mach
        super().__init__(mach)
        self.index = len(mach.lvars)
        ins.vars_index.append(self.index)
        mach.lvars.append(self)
        self.init_value = init_value
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def to_exp(self):
        e = Expression.id_lvar(self.mach, self.index)
        return e
