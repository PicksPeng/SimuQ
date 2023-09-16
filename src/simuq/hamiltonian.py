"""
The class file for time independent Hamiltonian.

A TIHamiltonian is stored as a sum of product terms:
a coefficient times a product of basic operators.
The products of basic operators are stored as a list
of operators (currently strings).

In the worst case, this representation may take
exponential resources. In practice, it is enough for
nowadays machines.

Natural computation operations are overloaded.
Pauli basis products are symbolically calculated, so
that commutativity test is possible.
"""

from copy import copy

from simuq.expression import Expression


class TIHamiltonian:
    """The time-independent Hamiltonian

    The underlying data-structure to store it is a list of
    tuples (h, c). Here h is a product Hamiltonian,
    represented by a list of site operators on the sites.
    c is the coefficient corresponding to h, and can either
    be a real number or an Expression.

    The site operator algebras are symbolically dealt with
    here, since fermionic operators' algebra is position-related.

    We can also test commutativity of Hamiltonians, and
    calculate the sites the Hamiltonian is acting non-
    trivially on.
    """

    def __init__(self, sites_type, ham):
        self.sites_type = sites_type
        self.ham = ham
        self.saved_t_sites = None
        self.cleanHam()

    @classmethod
    def empty(cls, sites_type):
        ham = []
        return cls(sites_type, ham)

    @classmethod
    def identity(cls, sites_type):
        num_sites = len(sites_type)
        prod = ["" for i in range(num_sites)]
        ham = [(prod, 1)]
        return cls(sites_type, ham)

    @classmethod
    def op(cls, sites_type, index, op):
        num_sites = len(sites_type)
        prod = ["" for i in range(num_sites)]
        prod[index] = op
        ham = [(prod, 1)]
        return cls(sites_type, ham)

    def operAlgebra(self):
        self.extend_ham_by_sites()

        i = 0
        while i < len(self.ham):
            (h, coef) = self.ham[i]
            for j in range(len(h)):
                if self.sites_type[j] == "qubit":
                    if h[j] == "XX" or h[j] == "YY" or h[j] == "ZZ":
                        h[j] = ""
                    elif h[j] == "XY":
                        h[j] = "Z"
                        coef *= 1j
                    elif h[j] == "YX":
                        h[j] = "Z"
                        coef *= -1j
                    elif h[j] == "YZ":
                        h[j] = "X"
                        coef *= 1j
                    elif h[j] == "ZY":
                        h[j] = "X"
                        coef *= -1j
                    elif h[j] == "ZX":
                        h[j] = "Y"
                        coef *= 1j
                    elif h[j] == "XZ":
                        h[j] = "Y"
                        coef *= -1j
            for j in range(len(h)):
                if self.sites_type[j] == "boson":
                    s = h[j].find("ac")
                    while s != -1:
                        h_copy = h.copy()
                        h_copy[j] = h_copy[j].replace("ac", "ca", 1)
                        h[j] = h[j].replace("ac", "", 1)
                        self.ham.append((h_copy, coef))
                        s = h[j].find("ac")
                elif self.sites_type[j] == "fermion":
                    s = h[j].find("ac")
                    while s != -1:
                        h_copy = h.copy()
                        h_copy[j] = h_copy[j].replace("ac", "ca", 1)
                        h[j] = h[j].replace("ac", "", 1)
                        self.ham.append((h_copy, -coef))
                        s = h[j].find("ac")
            self.ham[i] = (h, coef)
            i += 1

    def cleanHam(self, tol=1e-6):
        self.operAlgebra()

        hamdic = dict([])
        for h, coef in self.ham:
            htup = tuple(h)
            if htup not in hamdic:
                hamdic[htup] = coef
            else:
                hamdic[htup] += coef

        self.ham = []
        for htup in hamdic:
            if isinstance(hamdic[htup], Expression) or abs(hamdic[htup]) > tol:
                self.ham.append((list(htup), hamdic[htup]))

        """
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
            if abs(self.ham[i][1]) < tol :
                del self.ham[i]
                continue
            i += 1
        """

    def extend_ham_by_sites(self):
        for i in range(len(self.ham)):
            (h, coef) = self.ham[i]
            if len(h) < len(self.sites_type):
                h += [""] * (len(self.sites_type) - len(h))
                self.ham[i] = (h, coef)

    def extend_sites(self, new_sites_type):
        if len(new_sites_type) <= len(self.sites_type):
            return
        for i in range(len(self.sites_type)):
            if self.sites_type[i] != new_sites_type[i]:
                raise Exception("Sites type inconsistent!")
        add_num = len(new_sites_type) - len(self.sites_type)
        new_ham = [(self.ham[i][0] + [""] * add_num, self.ham[i][1]) for i in range(len(self.ham))]
        self.sites_type = new_sites_type
        self.ham = new_ham
        if self.saved_t_sites is not None:
            self.saved_t_sites += [0 for j in range(add_num)]

    def __neg__(self):
        ham = copy(self.ham)
        for i in range(len(ham)):
            ham[i] = (ham[i][0], -ham[i][1])
        h = TIHamiltonian(self.sites_type, ham)
        return h

    def __add__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            other = other * TIHamiltonian.identity(self.sites_type)
        if self.sites_type != other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        ham = copy(self.ham)
        ham += other.ham
        h = TIHamiltonian(self.sites_type, ham)
        h.cleanHam()
        return h

    def __radd__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return self.__add__(other * TIHamiltonian.identity(self.sites_type))
        else:
            return NotImplemented

    def __sub__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            other = other * TIHamiltonian.identity(self.sites_type)
        if self.sites_type != other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        ham = copy(self.ham)
        ham += other.__neg__().ham
        h = TIHamiltonian(self.sites_type, ham)
        h.cleanHam()
        return h

    def __rsub__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return (other * TIHamiltonian.identity(self.sites_type)) - self
        else:
            return NotImplemented

    @staticmethod
    def strlist_mul(sites_type, a, b):
        c = []
        coef = 1

        # Use a suffix sum to calculate the number of fermionic operators after site i.
        num_ferm_ope_backward = [0 for i in range(len(sites_type))]
        for i in range(len(sites_type) - 1, 0, -1):
            num_ferm_ope_backward[i - 1] = num_ferm_ope_backward[i]
            if sites_type[i] == "fermion":
                num_ferm_ope_backward[i - 1] += len(a[i])

        for i in range(len(sites_type)):
            if sites_type[i] == "qubit":
                c.append(a[i] + b[i])
            elif sites_type[i] == "boson":
                c.append(a[i] + b[i])
            elif sites_type[i] == "fermion":
                c.append(a[i] + b[i])
                if num_ferm_ope_backward[i] % 2 == 1 and len(b[i]) % 2 == 1:
                    coef *= -1

        return (c, coef)

    def scalar_mul(self, other):
        ham = copy(self.ham)
        for i in range(len(ham)):
            ham[i] = (ham[i][0], ham[i][1] * other)
        h = TIHamiltonian(self.sites_type, ham)
        return h

    def __mul__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            return self.scalar_mul(other)
        self.extend_ham_by_sites()
        other.extend_ham_by_sites()
        if self.sites_type != other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        ham = []
        for prod1, coef1 in self.ham:
            for prod2, coef2 in other.ham:
                (prod, coef3) = self.strlist_mul(self.sites_type, prod1, prod2)
                ham.append((prod, coef1 * coef2 * coef3))
        h = TIHamiltonian(self.sites_type, ham)
        h.cleanHam()
        return h

    def __truediv__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return self.scalar_mul(1 / other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return self.scalar_mul(other)
        else:
            return NotImplemented

    def exp_eval(self, gvars, lvars):
        ham = copy(self.ham)
        for i in range(len(ham)):
            ham[i] = (ham[i][0], ham[i][1].exp_eval(gvars, lvars))
        h = TIHamiltonian(self.sites_type, ham)
        return h

    def is_empty(self):
        abs_sum = 0
        for prod, c in self.ham:
            abs_sum += abs(c)
        return abs_sum < 1e-8

    @staticmethod
    def commutativity_test(h1, h2, derived=False):
        if derived:
            t_sites1 = h1.touched_sites()
            t_sites2 = h2.touched_sites()
            for i in range(len(h1.sites_type)):
                if t_sites1[i] == 1 and t_sites2[i] == 1:
                    return False
            return True
        return (h1 * h2 - h2 * h1).is_empty()

    def touched_sites(self):
        if self.saved_t_sites is not None:
            return self.saved_t_sites
        ret = [0 for i in range(len(self.sites_type))]
        for prod, t in self.ham:
            for i in range(len(self.sites_type)):
                if prod[i] != "":
                    ret[i] = 1
        self.saved_t_sites = ret
        return ret

    def to_qiskit_opflow(self):
        from qiskit.opflow import I, X, Y, Z

        def strlist_to_oplist(l):
            ret = []
            for i in range(len(l)):
                if l[i] == "":
                    ret.append(I)
                elif l[i] == "X":
                    ret.append(X)
                elif l[i] == "Y":
                    ret.append(Y)
                else:
                    ret.append(Z)
            return ret

        def list_kron(l):
            res = l[0]
            for i in range(1, len(l)):
                res = res ^ l[i]
            return res

        h = 0
        for prod, c in self.ham:
            h += c * list_kron(strlist_to_oplist(prod))

        return h

    def to_qutip_qobj(self):
        from qutip import qeye, sigmax, sigmay, sigmaz, tensor

        def strlist_to_oplist(l):
            ret = []
            for i in range(len(l)):
                if l[i] == "":
                    ret.append(qeye(2))
                elif l[i] == "X":
                    ret.append(sigmax())
                elif l[i] == "Y":
                    ret.append(sigmay())
                else:
                    ret.append(sigmaz())
            return ret

        ret = 0
        for prod, c in self.ham:
            ret = ret + tensor(strlist_to_oplist(prod)) * c

        return ret
