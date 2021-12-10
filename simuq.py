from copy import deepcopy
import scipy.optimize as opt
import numpy as np
import networkx as nx
import math

class BaseVar :
    def __init__(self, mach) :
        self.mach = mach

    def to_exp(self) :
        pass

    def __neg__(self) :
        return -self.to_exp()

    def __add__(self, other) :
        return self.to_exp() + other

    def __mul__(self, other) :
        return self.to_exp() * other

    def __sub__(self, other) :
        return self.to_exp() - other

    def __truediv__(self, other) :
        return self.to_exp() / other

class Expression :
    def __init__(self, mach, exp) :
        self.mach = mach
        self.exp = exp

    @classmethod
    def unit(cls, mach) :
        exp = lambda gvars, lvars : 1
        return cls(mach, exp)
    
    @classmethod
    def id_gvar(cls, mach, index) :
        exp = lambda gvars, lvars : gvars[index]
        return cls(mach, exp)

    @classmethod
    def id_lvar(cls, mach, index) :
        exp = lambda gvars, lvars : lvars[index]
        return cls(mach, exp)

    @classmethod
    def cos(cls, e) :
        if type(e) == int  or  type(e) == float  or  type(e) == complex :
            return math.cos(e)
        if isinstance(e, BaseVar) :
            e = e.to_exp()
        exp = lambda gvars, lvars : math.cos(e.exp(gvars, lvars))
        return cls(e.mach, exp)

    def __neg__(self) :
        exp = lambda gvars, lvars : -self.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def __add__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            exp = lambda gvars, lvars : self.exp(gvars, lvars) + other
            e = Expression(self.mach, exp)
            return e
        if isinstance(other, BaseVar) :
            other = other.to_exp()
        if not hasattr(other, "exp") :
            return NotImplemented
        exp = lambda gvars, lvars : self.exp(gvars, lvars) + other.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def __sub__(self, other) :
        return self.__add__(other.__neg__())

    def __radd__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            return self.__add__(other)
        elif isinstance(other, BaseVar) :
            other = other.to_exp()
            return self.__add__(other)
        else :
            return NotImplemented

    def __mul__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            exp = lambda gvars, lvars : self.exp(gvars, lvars) * other
            e = Expression(self.mach, exp)
            return e
        if isinstance(other, BaseVar) :
            other = other.to_exp()
        if not hasattr(other, "exp") :
            return NotImplemented
        exp = lambda gvars, lvars : self.exp(gvars, lvars) * other.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def __rmul__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            return self.__mul__(other)
        elif isinstance(other, BaseVar) :
            other = other.to_exp()
            return self.__mul__(other)
        else :
            return NotImplemented

    def __pow__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            exp = lambda gvars, lvars : self.exp(gvars, lvars) ** other
            e = Expression(self.mach, exp)
            return e
        if isinstance(other, BaseVar) :
            other = other.to_exp()
        if not hasattr(other, "exp") :
            return NotImplemented
        exp = lambda gvars, lvars : self.exp(gvars, lvars) ** other.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def __truediv__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            exp = lambda gvars, lvars : self.exp(gvars, lvars) / other
            e = Expression(self.mach, exp)
            return e
        if isinstance(other, BaseVar) :
            other = other.to_exp()
        if not hasattr(other, "exp") :
            return NotImplemented
        exp = lambda gvars, lvars : self.exp(gvars, lvars) / other.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def __rtruediv__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex :
            return self.__mul__(other)
        elif isinstance(other, BaseVar) :
            other = other.to_exp()
            return self.__mul__(other)
        else :
            return NotImplemented

    def exp_eval(self, gvars, lvars) :
        return self.exp(gvars, lvars)

class TIHamiltonian :
    def __init__(self, num_sites, ham) :
        self.num_sites = num_sites
        self.ham = ham
        self.saved_t_sites = None

    @classmethod
    def empty(cls, num_sites) :
        ham = []
        return cls(num_sites, ham)

    @classmethod
    def identity(cls, num_sites) :
        prod = ['' for i in range(num_sites)]
        ham = [(prod, 1)]
        return cls(num_sites, ham)

    @classmethod
    def op(cls, num_sites, index, op) :
        prod = ['' for i in range(num_sites)]
        prod[index] = op
        ham = [(prod, 1)]
        return cls(num_sites, ham)

    def cleanHam(self) :
        i = 0
        while i < len(self.ham) :
            for j in range(i) :
                if self.ham[j][0] == self.ham[i][0] :
                    self.ham[j] = (self.ham[j][0], self.ham[j][1] + self.ham[i][1])
                    del self.ham[i]
                    i -= 1
                    break
            i += 1
        i = 0
        while i < len(self.ham) :
            if isinstance(self.ham[i][1], Expression) :
                break
            if abs(self.ham[i][1]) < 1e-8 :
                del self.ham[i]
                continue
            i += 1

    def __neg__(self) :
        ham = deepcopy(self.ham)
        for i in range(len(ham)) :
            ham[i] = (ham[i][0], -ham[i][1])
        h = TIHamiltonian(self.num_sites, ham)
        return h

    def __add__(self, other) :
        # need more typing restrictions
        if self.num_sites != other.num_sites :
            return NotImplemented
        ham = deepcopy(self.ham)
        ham += other.ham
        h = TIHamiltonian(self.num_sites, ham)
        h.cleanHam()
        return h

    def __sub__(self, other) :
        if self.num_sites != other.num_sites :
            return NotImplemented
        ham = deepcopy(self.ham)
        ham += other.__neg__().ham
        h = TIHamiltonian(self.num_sites, ham)
        h.cleanHam()
        return h

    @staticmethod
    def strlist_mul(a, b) :
        if len(a) != len(b) :
            return NotImplementedError
        c = [a[i] + b[i] for i in range(len(a))]
        coef = 1
        for i in range(len(c)) :
            if c[i] == 'XX'  or  c[i] == 'YY'  or  c[i] == 'ZZ' :
                c[i] = ''
            elif c[i] == 'XY' :
                c[i] = 'Z'
                coef *= 1j
            elif c[i] == 'YX' :
                c[i] = 'Z'
                coef *= -1j
            elif c[i] == 'YZ' :
                c[i] = 'X'
                coef *= 1j
            elif c[i] == 'ZY' :
                c[i] = 'X'
                coef *= -1j
            elif c[i] == 'ZX' :
                c[i] = 'Y'
                coef *= 1j
            elif c[i] == 'XZ' :
                c[i] = 'Y'
                coef *= -1j
        return (c, coef)

    def scalar_mul(self, other) :
        ham = deepcopy(self.ham)
        for i in range(len(ham)) :
            ham[i] = (ham[i][0], ham[i][1] * other)
        h = TIHamiltonian(self.num_sites, ham)
        return h

    def __mul__(self, other) :
        # need more typing restrictions
        if type(other) == int  or  type(other) == float  or  type(other) == complex  or  isinstance(other, Expression) :
            return self.scalar_mul(other)
        if self.num_sites != other.num_sites :
            return NotImplemented
        ham = []
        for (prod1, coef1) in self.ham :
            for (prod2, coef2) in other.ham :
                (prod, coef3) = self.strlist_mul(prod1, prod2)
                ham.append((prod, coef1 * coef2 * coef3))
        h = TIHamiltonian(self.num_sites, ham)
        h.cleanHam()
        return h

    def __truediv__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex  or  isinstance(other, Expression) :
            return self.scalar_mul(1/other)
        else :
            return NotImplemented

    def __rmul__(self, other) :
        if type(other) == int  or  type(other) == float  or  type(other) == complex  or  isinstance(other, Expression) :
            return self.scalar_mul(other)
        else :
            return NotImplemented

    def exp_eval(self, gvars, lvars) :
        ham = deepcopy(self.ham)
        for i in range(len(ham)) :
            ham[i] = (ham[i][0], ham[i][1].exp_eval(gvars, lvars))
        h = TIHamiltonian(self.num_sites, ham)
        return h

    def is_empty(self) :
        abs_sum = 0
        for (prod, c) in self.ham :
            abs_sum += abs(c)
        if abs_sum < 1e-8 :
            return True
        else :
            return False

    @staticmethod
    def commutativity_test(h1, h2, derived = False) :
        if derived :
            t_sites1 = h1.touched_sites()
            t_sites2 = h2.touched_sites()
            for i in range(h1.num_sites) :
                if t_sites1[i] == 1  and  t_sites2[i] == 1 :
                    return False
            return True
        return (h1 * h2 - h2 * h1).is_empty()

    def touched_sites(self) :
        if self.saved_t_sites != None :
            return self.saved_t_sites
        ret = [0 for i in range(self.num_sites)]
        for (prod, t) in self.ham :
            for i in range(self.num_sites) :
                if prod[i] != '' :
                    ret[i] = 1
        self.saved_t_sites = ret
        return ret


class BaseQuantumEnvironment :
    def __init__(self) :
        self.num_sites = 0
        self.sites = []
        self.sites_type = []

    def identity(self) :
        return TIHamiltonian.identity(self.num_sites)

    def singletonOp(self, index, op) :
        return TIHamiltonian.op(self.num_sites, index, op)




class QSystem(BaseQuantumEnvironment) :
    def __init__(self) :
        super().__init__()
        self.evos = []

    def addEvolution(self, h, t) :
        self.evos.append((h, t))



class BaseSite :
    def __init__(self, qs) :
        self.index = qs.num_sites
        self.qs = qs
        qs.num_sites += 1
        qs.sites.append(self)

    def createOp(self, op) :
        h = self.qs.singletonOp(self.index, op)
        return h
        

class qubit(BaseSite) :
    def __init__(self, qs) :
        super().__init__(qs)
        qs.sites_type.append('qubit')

    def X(self) :
        return self.createOp("X")

    def Y(self) :
        return self.createOp("Y")

    def Z(self) :
        return self.createOp("Z")

class fock(BaseSite) :
    def __init__(self, qs) :
        super().__init__(qs)
        qs.sites_type.append('fock')

    def a(self) :
        return self.createOp("a")

    def c(self) :
        return self.createOp("c")


# Jaynes-Cummings model
qs = QSystem()
atom = qubit(qs)
cav = fock(qs)

omega_a = 0.1
omega_c = 0.01
g = 0.001
h = omega_a * atom.Z() / 2 + omega_c * cav.c() * cav.a() + g * (cav.a() + cav.c()) * atom.X()
t = 0.1
qs.addEvolution(h, t)











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
        newh = deepcopy(h)
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



Rydberg = QMachine()

q1 = qubit(Rydberg)
q2 = qubit(Rydberg)
q3 = qubit(Rydberg)

x0 = 0
x1 = GlobalVar(Rydberg)
x2 = GlobalVar(Rydberg)

L1 = SignalLine(Rydberg)

ins1 = Instruction(L1, 'native', 'ins1')
a = LocalVar(ins1)
ins1.set_ham((Expression.cos(a) + 2 / (x2 - x1) ** 4) * q1.X() * q2.X())

ins2 = Instruction(L1, 'derived', 'ins2')
a = LocalVar(ins2)
ins2.set_ham((a * (x1 + x0)) * (q2.X() + q1.Y()))





def locate_switch_evo(mach, evo_index) :
    return mach.num_gvars + evo_index * (mach.num_inss + mach.num_lvars)

def switch_term(mach, evo_index, ins_index, mc) :
    return (lambda x : x[locate_switch_evo(mach, evo_index) + ins_index] * mc.exp_eval(x[:mach.num_gvars], x[locate_switch_evo(mach, evo_index) + mach.num_inss : locate_switch_evo(mach, evo_index + 1)]))

def switch_fun(mach, evo_index, ins_index) :
    return (lambda x : x[locate_switch_evo(mach, evo_index) + ins_index])

def locate_nonswitch_evo(mach, evo_index) :
    return mach.num_gvars + evo_index * mach.num_lvars

def non_switch_term(mach, evo_index, ins_index, mc) :
    return (lambda x : mc.exp_eval(x[:mach.num_gvars], x[locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(mach, evo_index + 1)]))

def solve_aligned(ali, qs, mach) :
    print(ali)
    
    eqs = []
    nvar = mach.num_gvars + len(qs.evos) * (mach.num_inss + mach.num_lvars)
    for evo_index in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_index]
        mark = [[0 for ins in line.inss] for line in mach.lines]
        for (tprod, tc) in h.ham :
            eq = (lambda c : lambda x : -c)(tc)
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    for (mprod, mc) in ins.h.ham :
                        prod_match = True
                        for k in range(qs.num_sites) :
                            if tprod[k] != mprod[ali[k]] :
                                prod_match = False
                                break
                        if prod_match :
                            eq = (lambda eq_, f_ : lambda x : eq_(x) + f_(x))(eq, switch_term(mach, evo_index, ins.index, mc))
                            mark[i][j] = 1
                            break
            eqs.append(eq)
        for i in range(len(mach.lines)) :
            line = mach.lines[i]
            for j in range(len(line.inss)) :
                ins = line.inss[j]
                if mark[i][j] == 0 :
                    eqs.append(switch_fun(mach, evo_index, ins.index))

    #if len(eqs) < nvar :
    #    eqs = eqs + [(lambda x : 0) for i in range(nvar - len(eqs))]
    f = (lambda eqs_ : lambda x : [(lambda i_ : eqs_[i_](x))(i) for i in range(len(eqs_))])(eqs)
    var_lb = -np.inf
    var_ub = np.inf
    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    init = [0 for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [0 for j in range(mach.num_inss)] + [var_lb for j in range(mach.num_lvars)]
        ubs += [1 for j in range(mach.num_inss)] + [var_ub for j in range(mach.num_lvars)]
        init += [0.5 for j in range(mach.num_inss)] + [0 for j in range(mach.num_lvars)]

    sol_detail = opt.least_squares(f, init, bounds = (lbs, ubs))
    sol = sol_detail.x

    print(np.linalg.norm(f(sol)))
    if np.linalg.norm(f(sol)) > 1e-3 :
        return False

    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [var_lb for j in range(mach.num_lvars)]
        ubs += [var_ub for j in range(mach.num_lvars)]

    eqs = []
    nvar = mach.num_gvars + len(qs.evos) * mach.num_lvars
    initvars = sol[:mach.num_gvars].tolist()
    switch = [[[0 for ins in line.inss] for line in mach.lines] for evo_index in range(len(qs.evos))]
    for evo_index in range(len(qs.evos)) :
        for i in range(len(mach.lines)) :
            line = mach.lines[i]
            for j in range(len(line.inss)) :
                ins = line.inss[j]
                if abs(sol[locate_switch_evo(mach, evo_index) + ins.index]) < 1e-3 :
                    switch[evo_index][i][j] = 0
                else :
                    switch[evo_index][i][j] = 1
        initvars += sol[locate_switch_evo(mach, evo_index) + mach.num_inss : locate_switch_evo(mach, evo_index + 1)].tolist()
    for evo_index in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_index]
        for (tprod, tc) in h.ham :
            eq = (lambda c : lambda x : -c)(tc)
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    if switch[evo_index][i][j] == 1 :
                        for (mprod, mc) in ins.h.ham :
                            prod_match = True
                            for k in range(qs.num_sites) :
                                if tprod[k] != mprod[ali[k]] :
                                    prod_match = False
                                    break
                            if prod_match :
                                eq = (lambda eq_, f_ : lambda x : eq_(x) + f_(x))(eq, non_switch_term(mach, evo_index, ins.index, mc))
                                break
            eqs.append(eq)

    f = (lambda eqs_ : lambda x : [(lambda i_ : eqs_[i_](x))(i) for i in range(len(eqs_))])(eqs)
    var_lb = -np.inf
    var_ub = np.inf
    lbs = [var_lb for i in range(mach.num_gvars)]
    ubs = [var_ub for i in range(mach.num_gvars)]
    for i in range(len(qs.evos)) :
        lbs += [var_lb for j in range(mach.num_lvars)]
        ubs += [var_ub for j in range(mach.num_lvars)]

    sol_detail = opt.least_squares(f, initvars, bounds = (lbs, ubs))
    sol = sol_detail.x

    global gsol
    gsol = sol
    global gswitch
    gswitch = switch
    
    print(np.linalg.norm(f(sol)))
    if np.linalg.norm(f(sol)) > 1e-5 :
        return False

    print(switch)
    print(sol)
    return True


def align(i, ali, qs, mach) :
    if i == qs.num_sites :
        if solve_aligned(ali, qs, mach) :
            return True
        return False
    for x in range(mach.num_sites) :
        available = True
        ali[i] = x
        for j in range(i) :
            if ali[j] == x :
                available = False
                break
        if available == False :
            continue
        for (h, t) in qs.evos :
            for (tprod, tc) in h.ham :
                found = False
                for line in mach.lines :
                    for ins in line.inss :
                        for (mprod, mc) in ins.h.ham :
                            prod_partial_match = True
                            for k in range(i + 1) :
                                if tprod[k] != mprod[ali[k]] :
                                    prod_partial_match = False
                            if prod_partial_match :
                                found = True
                                break
                        if found :
                            break
                    if found :
                        break
                if not found :
                    available = False
                    break
            if not available :
                break
        if available :
            if align(i + 1, ali, qs, mach) :
                return True
    return False


def find_sol(qs, mach) :
    return align(0, [0 for i in range(qs.num_sites)], qs, mach)


qs = QSystem()
q1 = qubit(qs)
q2 = qubit(qs)
c = fock(qs)
h = 0.5 * q2.X() * c.a() + 2 * q1.X() * c.c()
qs.addEvolution(h, 1)
h = 0.2 * q2.X() * c.a() + c.a() * c.c()
qs.addEvolution(h, 1)


mach = QMachine()
c = fock(mach)
q1 = qubit(mach)
q2 = qubit(mach)
q3 = qubit(mach)

L1 = SignalLine(mach)

ins1 = Instruction(L1, 'native', 'L1_ins1')
a = LocalVar(ins1)
ins1.set_ham(Expression.cos(a) * c.a() * q2.X())

ins2 = Instruction(L1, 'derived', 'L1_ins2')
a = LocalVar(ins2)
ins2.set_ham(a * q1.X() * c.c())

L2 = SignalLine(mach)

ins3 = Instruction(L2, 'native', 'L2_ins1')
ins3.set_ham(c.a() * c.c())


if find_sol(qs, mach) :
    trotter_step = 4
    sol = gsol
    switch = gswitch
    sol_gvars = sol[:mach.num_gvars]
    boxes = []
    edges = []
    ending_boxes = []
    for evo_index in range(len(qs.evos)) :
        (h, t) = qs.evos[evo_index]
        next_ending_boxes = []
        coloring = [i for i in range(mach.num_sites)]
        sol_lvars = sol[locate_nonswitch_evo(mach, evo_index) : locate_nonswitch_evo(mach, evo_index + 1)].tolist()
        # Detach by touched sites
        first_touch = [[-1 for i in range(len(line.inss))] for line in mach.lines]
        for i in range(len(mach.lines)) :
            line = mach.lines[i]
            for j in range(len(line.inss)) :
                ins = line.inss[j]
                if switch[evo_index][i][j] == 1 :
                    t_sites = (ins.h.exp_eval(sol_gvars, sol_lvars)).touched_sites()
                    for k in range(mach.num_sites) :
                        if t_sites[k] == 1 :
                            if first_touch[i][j] == -1 :
                                first_touch[i][j] = k
                            else :
                                if coloring[k] != coloring[first_touch[i][j]] :
                                    c = coloring[k]
                                    for l in range(mach.num_sites) :
                                        if coloring[l] == c :
                                            coloring[l] = coloring[first_touch[i][j]]
        color_part = []
        for k in range(mach.num_sites) :
            ins_set = []
            for i in range(len(mach.lines)) :
                line = mach.lines[i]
                for j in range(len(line.inss)) :
                    ins = line.inss[j]
                    if switch[evo_index][i][j] == 1  and  coloring[first_touch[i][j]] == k :
                        ins_lvars = []
                        for i in range(len(ins.vars_index)) :
                            ins_lvars.append(sol_lvars[ins.vars_index[i]])
                        ins_set.append(((i, j), ins, ins.exp_eval(sol_gvars, sol_lvars), ins_lvars))
            if ins_set != [] :
                color_part.append(ins_set)
        # Detach if commute with all others
        for i in range(len(color_part)) :
            ins_set = color_part[i]
            j = 0
            while j < len(ins_set) :
                all_commute = True
                for k in range(len(ins_set)) :
                    if j == k :
                        continue
                    if not TIHamiltonian.commutativity_test(ins_set[j][2], ins_set[k][2], \
                                                            ins_set[j][1].prop == 'derived' or ins_set[k][1].prop == 'derived') :
                        all_commute = False
                        break
                if all_commute :
                    color_part.append([ins_set[j]])
                    del ins_set[j]
                    j -= 1
                j += 1
        # Trotterization
        for ins_set in color_part :
            local_ending_boxes = ending_boxes
            n = len(ins_set)
            G = nx.Graph()
            G.add_nodes_from(range(n))
            for i in range(n) :
                for j in range(n) :
                    if i == j :
                        continue
                    if (ins_set[i][0][0] == ins_set[j][0][0]) \
                       or ((ins_set[i][1].prop == 'derived' or ins_set[j][1].prop == 'derived') \
                           and (not TIHamiltonian.commutativity_test(ins_set[i][2], ins_set[j][2], True))) :
                        G.add_edge(i, j)
            col = nx.coloring.greedy_color(G, strategy = 'largest_first')
            nodes_of_color = [[] for i in set(col.values())]
            for i in range(n) :
                nodes_of_color[col[i]].append(i)
            # Check if the partitions are pair-wise commute
            sumh = []
            for i in range(len(nodes_of_color)) :
                h = TIHamiltonian.empty(mach.num_sites)
                for j in range(len(nodes_of_color[i])) :
                    h = h + ins_set[nodes_of_color[i][j]][2]
                sumh.append(h)
            all_commute = True
            for i in range(len(nodes_of_color)) :
                for j in range(i) :
                    if TIHamiltonian.commutativity_test(sumh[i], sumh[j]) == False :
                        all_commute = False
                        break
                if not all_commute :
                    break
            if all_commute :
                steps = 1
            else :
                steps = trotter_step
            for i in range(steps) :
                next_local_ending_boxes = []
                for k in range(len(nodes_of_color)) :
                    list_ins = nodes_of_color[k]
                    box_label = len(boxes)
                    box_ins = []
                    for j in range(len(list_ins)) :
                        box_ins.append(ins_set[list_ins[j]])
                    boxes.append((box_ins, t / trotter_step, sumh[k]))
                    for label in local_ending_boxes :
                        edges.append((label, box_label))
                    next_local_ending_boxes.append(box_label)
                local_ending_boxes = next_local_ending_boxes
            next_ending_boxes += next_local_ending_boxes

        ending_boxes = next_ending_boxes

    # Delete commutative edges
    G = nx.DiGraph()
    G.add_nodes_from(range(len(boxes)))
    for edge in edges :
        if edge[1] not in nx.ancestors(G, edge[0]) :
            G.add_edge(*edge)
    s = 0
    while s < len(edges) :
        edge = edges[s]
        h1 = boxes[edge[0]][2]
        h2 = boxes[edge[1]][2]
        if TIHamiltonian.commutativity_test(h1, h2) :
            del edges[s]
            G.remove_edge(*edge)
            for i in range(len(edges)) :
                if edges[i][1] == edge[0] :
                    new_edge = (edges[i][0], edge[1])
                    if new_edge[1] not in nx.ancestors(G, new_edge[0]) :
                        edges.append(new_edge)
                        G.add_edge(*new_edge)
            s -= 1
        s += 1
    
    print()
    print(sol_gvars)
    for box in boxes :
        print(([((i, j), ins_lvars) for ((i, j), ins, h, ins_lvars) in box[0]], box[1]))
    for edge in edges :
        print(edge)
