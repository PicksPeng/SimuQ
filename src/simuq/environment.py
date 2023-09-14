"""
The classes describing quantum environments.

A quantum environment is a container of sites. Either target quantum
systems or quantum machines are described by its children classes.

A site is the basic unit describing a quantized physics entity. It
includes qubit / bosonic / fermionic states, or other customized
types of sites. Each site contains a set of site operators. These
operators will be used in the construction of Hamiltonians.

Currently operators are stored as strings. In future these may be
substituted by operator classes.
"""

from simuq.hamiltonian import TIHamiltonian


class BaseQuantumEnvironment:
    """The basic quantum environment.

    It models a system to which quantum sites belong to.

    The sites are stored in a sequence, where their types are stored in
    list sites_type.
    """

    def __init__(self):
        self.sites = []
        self.sites_type = []
        self.num_sites = 0

    def identity(self):
        return TIHamiltonian.identity(self.sites_type)

    def singletonOp(self, index, op):
        return TIHamiltonian.op(self.sites_type, index, op)


class BaseSite:
    """The basic quantum site."""

    def __init__(self, qs, name=None):
        self.index = qs.num_sites
        self.qs = qs
        if name is None:
            name = f"Site{qs.num_sites}"
        self.name = name
        qs.num_sites += 1
        qs.sites.append(self)

    def createOp(self, op):
        h = self.qs.singletonOp(self.index, op)
        return h


class Qubit(BaseSite):
    """The qubit site.

    By default, there are X, Y, Z, I defined as site operators
    of a qubit site.
    """

    def __init__(self, qs, name=None):
        if name is None:
            name = f"Qubit{qs.num_sites}"
        super().__init__(qs, name)
        qs.sites_type.append("qubit")
        self.X = self.gen_X()
        self.Y = self.gen_Y()
        self.Z = self.gen_Z()
        self.I = self.gen_I()

    def gen_X(self):
        return self.createOp("X")

    def gen_Y(self):
        return self.createOp("Y")

    def gen_Z(self):
        return self.createOp("Z")

    def gen_I(self):
        return self.createOp("")


class BaseParticle(BaseSite):
    """The basic particle site.

    By default, there are annihilation and creation operators.
    Additionally, to be consistent with qubit, I represents the identity.
    """

    def __init__(self, qs, name=None):
        if name is None:
            name = f"Site{qs.num_sites}"
        super().__init__(qs, name)
        qs.sites_type.append("particle")
        self.a = self.gen_a()
        self.c = self.gen_c()
        self.I = self.gen_I()

    def gen_a(self):
        return self.createOp("a")

    def gen_c(self):
        return self.createOp("c")

    def gen_I(self):
        return self.createOp("")


class Fermion(BaseParticle):
    """The fermionic site"""

    def __init__(self, qs, name=None):
        if name is None:
            name = f"Fermion{qs.num_sites}"
        super().__init__(qs, name)
        qs.sites_type[-1] = "fermion"


class Boson(BaseParticle):
    """The bosonic site"""

    def __init__(self, qs, name=None):
        if name is None:
            name = f"Boson{qs.num_sites}"
        super().__init__(qs, name)
        qs.sites_type[-1] = "boson"
