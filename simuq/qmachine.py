from copy import deepcopy, copy
from simuq.environment import BaseQuantumEnvironment
from simuq.expression import BaseVar, Expression

class QMachine(BaseQuantumEnvironment) :
    def __init__(self) :
        super().__init__()
        self.num_gvars = 0
        self.gvars = []
        self.num_lvars = 0
        self.lvars = []
        self.num_lines = 0
        self.lines = []
        self.num_inss = 0


class SignalLine :
    def __init__(self, mach) :
        self.mach = mach
        self.index = mach.num_lines
        mach.num_lines += 1
        mach.lines.append(self)
        self.inss = []


class Instruction :
    def __init__(self, line, prop, name = 'ins') :
        self.mach = line.mach
        self.line = line
        line.inss.append(self)
        self.index = self.mach.num_inss
        self.mach.num_inss += 1
        self.vars_index = []
        self.prop = prop
        self.name = name

    def set_ham(self, h) :
        newh = copy(h)
        for i in range(len(h.ham)) :
            if not isinstance(h.ham[i][1], Expression) :
                newh.ham[i] = (h.ham[i][0], h.ham[i][1] * Expression.unit(self.mach))
        self.h = newh

    def exp_eval(self, gvars, lvars) :
        return self.h.exp_eval(gvars, lvars)


class GlobalVar(BaseVar) :
    def __init__(self, mach) :
        super().__init__(mach)
        self.index = mach.num_gvars
        mach.num_gvars += 1
        mach.gvars.append(self)

    def to_exp(self) :
        e = Expression.id_gvar(self.mach, self.index)
        return e
        

class LocalVar(BaseVar) :
    def __init__(self, ins) :
        mach = ins.mach
        super().__init__(mach)
        self.index = mach.num_lvars
        ins.vars_index.append(self.index)
        mach.num_lvars += 1
        mach.lvars.append(self)

    def to_exp(self) :
        e = Expression.id_lvar(self.mach, self.index)
        return e
