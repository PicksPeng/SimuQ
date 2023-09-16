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

import itertools
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

    def instantiate(self):
        if not self.instantiated:
            if self.with_sys_ham:
                sys_line = self.add_signal_line()
                sys_ins = sys_line.add_instruction("native", "System Hamiltonian")
                sys_ins.set_ham(self.sys_ham)
                sys_ins.is_sys_ham = True

            # Add indices to global variables
            for index, gvar in enumerate(self.gvars):
                gvar.set_index(index)
            # Populate self.lvars and set their indices
            self.lvars = [lvar for line in self.lines for ins in line.inss for lvar in ins.lvars]
            for index, lvar in enumerate(self.lvars):
                lvar.set_index(len(self.gvars) + index)
            for line in self.lines:
                for ins in line.inss:
                    ins.set_vars_index([lvar.index - len(self.gvars) for lvar in ins.lvars])

            # Set indices for instructions
            for index, ins in enumerate(itertools.chain(*[line.inss for line in self.lines])):
                ins.set_index(index)
            self.num_inss = sum(len(line.inss) for line in self.lines)

            self.instantiated = True

    def extend_instruction_sites(self):
        for line in self.lines:
            for ins in line.inss:
                ins.h.extend_sites(self.sites_type)

    def add_signal_line(self):
        line = SignalLine()
        self.lines.append(line)
        return line

    def add_global_variable(self, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        var = BaseVar(init_value, lower_bound, upper_bound)
        self.gvars.append(var)
        return var


class SignalLine:
    """A signal line.

    It contains all instructions belonging to it.
    """

    def __init__(self):
        self.inss = []

    def add_instruction(self, nativeness="native", name="ins"):
        ins = Instruction(nativeness, name)
        self.inss.append(ins)
        return ins


class Instruction:
    """An instruction.

    It contains the local variables belonging to it, its
    property, and its instruction Hamiltonian.
    """

    def __init__(self, nativeness="native", name="ins"):
        self.index = None
        self.lvars = []
        # Stores the index of the local variables in the QMachine (not lvar.index itself)
        self.vars_index = None
        self.nativeness = nativeness
        self.name = name
        self.is_sys_ham = False
        self.h = None

    def set_index(self, index):
        self.index = index

    def set_vars_index(self, vars_index):
        self.vars_index = vars_index

    def set_ham(self, h):
        newh = copy(h)
        for i in range(len(h.ham)):
            if not isinstance(h.ham[i][1], Expression):
                newh.ham[i] = (h.ham[i][0], h.ham[i][1] * Expression.unit())
        self.h = newh

    def exp_eval(self, gvars, lvars):
        return self.h.exp_eval(gvars, lvars)

    def add_local_variable(self, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        var = BaseVar(init_value, lower_bound, upper_bound)
        self.lvars.append(var)
        return var
