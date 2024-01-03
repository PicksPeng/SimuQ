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

from copy import copy, deepcopy
from collections import MutableMapping

from simuq.expression import Expression

class productHamiltonian(MutableMapping):
    """ A defaultdict-like object to speed up local Hamiltonian calculation.
        Does not expand with queries.
    """
    def __init__(self, from_list = None) :
        self.default = ''
        self.d = dict()
        if from_list :
            for (k, v) in from_list :
                self.d[k] = v

    def __getitem__(self, key) :
        if key in self.d :
            return self.d[key]
        else :
            return self.default

    def __setitem__(self, key, newvalue):
        self.d[key] = newvalue
        if newvalue == self.default :
            del self.d[key]

    def __delitem__(self, key):
        del self.d[key]

    def __iter__(self):
        return iter(self.d)

    def __len__(self):
        return len(self.d)

    def to_list(self) :
        keys = list(self.d.keys())
        keys.sort()
        l = []
        for k in keys :
            l.append((k, self.d[k]))
        return l

    def __eq__(self, other) :
        keys = set(self.d.keys()).union(other.d.keys())
        for k in keys :
            if self.d[k] != other.d[k] :
                return False
        return True

    def keys(self) :
        return self.d.keys()

    def __repr__(self) :
        return f"<prodHam d={self.d}>"


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
        #self.cleanHam()

    @classmethod
    def empty(cls, sites_type):
        ham = []
        return cls(sites_type, ham)

    @classmethod
    def identity(cls, sites_type):
        #num_sites = len(sites_type)
        #prod = ["" for i in range(num_sites)]
        prod = productHamiltonian()
        ham = [(prod, 1)]
        return cls(sites_type, ham)

    @classmethod
    def op(cls, sites_type, index, op):
        #num_sites = len(sites_type)
        #prod = ["" for i in range(num_sites)]
        prod = productHamiltonian()
        prod[index] = op
        ham = [(prod, 1)]
        return cls(sites_type, ham)

    def operAlgebra(self):
        #self.extend_ham_by_sites()

        i = 0
        while i < len(self.ham):
            (h, coef) = self.ham[i]
            #for j in range(len(h)):
            keys = list(h.keys())
            for j in keys:
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
            #for j in range(len(h)):
            keys = list(h.keys())
            for j in keys:
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
            #htup = tuple(h)
            htup = tuple(h.to_list())
            if htup not in hamdic:
                hamdic[htup] = coef
            else:
                hamdic[htup] += coef

        self.ham = []
        for htup in hamdic:
            if isinstance(hamdic[htup], Expression) or abs(hamdic[htup]) > tol:
                #self.ham.append((list(htup), hamdic[htup]))
                self.ham.append((productHamiltonian(from_list=htup), hamdic[htup]))

    """
    def extend_ham_by_sites(self):
        for i in range(len(self.ham)):
            (h, coef) = self.ham[i]
            if len(h) < len(self.sites_type):
                h += [""] * (len(self.sites_type) - len(h))
                self.ham[i] = (h, coef)
    """

    def extend_sites(self, new_sites_type):
        if new_sites_type is self.sites_type :
            return
        if new_sites_type == self.sites_type :
            return
        if len(new_sites_type) <= len(self.sites_type):
            return
        """ Delete type check for efficiency
        for i in range(len(self.sites_type)):
            if self.sites_type[i] != new_sites_type[i]:
                raise Exception("Sites type inconsistent!")
        """
        add_num = len(new_sites_type) - len(self.sites_type)
        self.sites_type = new_sites_type
        """
        new_ham = [(self.ham[i][0] + [""] * add_num, self.ham[i][1]) for i in range(len(self.ham))]
        self.ham = new_ham
        """
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
        if self.sites_type is not other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        ham = copy(self.ham)
        ham += other.ham
        h = TIHamiltonian(self.sites_type, ham)
        #h.cleanHam()
        return h

    def __radd__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return self.__add__(other * TIHamiltonian.identity(self.sites_type))
        else:
            return NotImplemented

    def __iadd__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            other = other * TIHamiltonian.identity(self.sites_type)
        if self.sites_type != other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        self.ham += other.ham
        return self

    def __sub__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            other = other * TIHamiltonian.identity(self.sites_type)
        if self.sites_type is not other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        ham = copy(self.ham)
        ham += other.__neg__().ham
        h = TIHamiltonian(self.sites_type, ham)
        #h.cleanHam()
        return h

    def __rsub__(self, other):
        if isinstance(other, (int, float, complex, Expression)):
            return (other * TIHamiltonian.identity(self.sites_type)) - self
        else:
            return NotImplemented

    def __isub__(self, other):
        # need more typing restrictions
        if isinstance(other, (int, float, complex, Expression)):
            other = other * TIHamiltonian.identity(self.sites_type)
        if self.sites_type != other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        self.ham += other.__neg__().ham
        return self

    @staticmethod
    def strlist_mul(sites_type, a, b):

        """
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
        """

        keys = list(set(a.keys()).union(b.keys()))
        keys.sort()
        num_ferm_ope_backward = [0 for i in range(len(keys))]
        for i in range(len(keys) - 1, 0, -1):
            num_ferm_ope_backward[i - 1] = num_ferm_ope_backward[i]
            if sites_type[keys[i]] == "fermion":
                num_ferm_ope_backward[i - 1] += len(a[keys[i]])

        c = productHamiltonian()
        coef = 1
        
        for i in keys:
            if sites_type[i] == "qubit":
                c[i] = a[i] + b[i]
            elif sites_type[i] == "boson":
                c[i] = a[i] + b[i]
            elif sites_type[i] == "fermion":
                c[i] = a[i] + b[i]
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
        #self.extend_ham_by_sites()
        #other.extend_ham_by_sites()
        if self.sites_type is not other.sites_type:
            self.extend_sites(other.sites_type)
            other.extend_sites(self.sites_type)
        self.cleanHam()
        other.cleanHam()
        ham = []
        for prod1, coef1 in self.ham:
            for prod2, coef2 in other.ham:
                (prod, coef3) = self.strlist_mul(self.sites_type, prod1, prod2)
                ham.append((prod, coef1 * coef2 * coef3))
        h = TIHamiltonian(self.sites_type, ham)
        h.cleanHam()
        return h

    def __pow__(self, exponent):
        if exponent < 0:
            raise ValueError("Negative exponent is not supported")
        if exponent == 0:
            return 1
        result = self
        for _ in range(exponent - 1):
            result = result * self
        return result

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
            """
            for i in range(len(self.sites_type)):
                if prod[i] != "":
                    ret[i] = 1
            """
            for i in prod:
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
            h += c * list_kron(strlist_to_oplist(prod.to_list()))

        return h

    def to_qutip_qobj(self):
        from qutip import qeye, sigmax, sigmay, sigmaz, tensor

        def strlist_to_oplist(l, n):
            ret = []
            for i in range(n):
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
            ret = ret + tensor(strlist_to_oplist(prod, len(self.sites_type))) * c

        return ret

    # A faster version of sum a list of TIHamiltonian
    @staticmethod
    def hlist_sum(hlist):
        n = len(hlist)
        if n == 0:
            return 0
        sites_type = hlist[0].sites_type
        for i in range(n):
            if hlist[i].sites_type != sites_type:
                raise Exception("Site types do not match")
        ham = []
        for i in range(n):
            ham += hlist[i].ham

        return TIHamiltonian(sites_type, ham)

hlist_sum = TIHamiltonian.hlist_sum
