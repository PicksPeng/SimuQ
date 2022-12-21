'''
The classes describing quantum environments.

A quantum environment is a container of sites. Either target quantum
systems or quantum machines are described by its children classes.

A site is the basic unit describing a quantized physics entity. It
includes qubit / fork states, or other customized types like qutrit
or qudit. One should provide sites with their supported operators.
These operators will be used in the construction of Hamiltonians.

Currently operators are stored as strings. In future these may be
substituted by operator classes.
'''

from simuq.hamiltonian import TIHamiltonian

class BaseQuantumEnvironment :
    def __init__(self) :
        self.sites = []
        self.sites_type = []
        self.num_sites = 0

    def identity(self) :
        return TIHamiltonian.identity(self.sites_type)

    def singletonOp(self, index, op) :
        return TIHamiltonian.op(self.sites_type, index, op)

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
        self.X = self.gen_X()
        self.Y = self.gen_Y()
        self.Z = self.gen_Z()
        self.I = self.gen_I()

    def gen_X(self) :
        return self.createOp("X")

    def gen_Y(self) :
        return self.createOp("Y")

    def gen_Z(self) :
        return self.createOp("Z")

    def gen_I(self) :
        return self.createOp("")


class BaseParticle(BaseSite) :
    def __init__(self, qs) :
        super().__init__(qs)
        qs.sites_type.append('particle')
        self.a = self.gen_a()
        self.c = self.gen_c()
        self.I = self.gen_I()

    def gen_a(self) :
        return self.createOp("a")

    def gen_c(self) :
        return self.createOp("c")

    def gen_I(self) :
        return self.createOp("")


class fermion(BaseParticle) :
    def __init__(self, qs) :
        super().__init__(qs)
        qs.sites_type[-1] = 'fermion'


class boson(BaseParticle) :
    def __init__(self, qs) :
        super().__init__(qs)
        qs.sites_type[-1] = 'boson'
    
