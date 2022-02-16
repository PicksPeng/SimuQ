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
        self.num_sites = 0
        self.sites = []
        self.sites_type = []

    def identity(self) :
        return TIHamiltonian.identity(self.num_sites)

    def singletonOp(self, index, op) :
        return TIHamiltonian.op(self.num_sites, index, op)

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
