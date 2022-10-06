'''
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
'''

from copy import copy, deepcopy
from simuq.expression import Expression

class TIHamiltonian :
    def __init__(self, num_sites, ham) :
        self.num_sites = num_sites
        self.ham = ham
        self.saved_t_sites = None

    @classmethod
    def empty(cls) :
        ham = []
        return cls(0, ham)

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

    def extend_sites(self, new_num_sites) :
        if new_num_sites <= self.num_sites :
            return
        add_num = new_num_sites - self.num_sites
        new_ham = [(self.ham[i][0] + ['' for j in range(add_num)], self.ham[i][1]) for i in range(len(self.ham))]
        self.num_sites = new_num_sites
        self.ham = new_ham
        if self.saved_t_sites is not None :
            self.saved_t_sites += [0 for j in range(add_num)]            

    def __neg__(self) :
        ham = copy(self.ham)
        for i in range(len(ham)) :
            ham[i] = (ham[i][0], -ham[i][1])
        h = TIHamiltonian(self.num_sites, ham)
        return h

    def __add__(self, other) :
        # need more typing restrictions
        if self.num_sites != other.num_sites :
            self.extend_sites(other.num_sites)
            other.extend_sites(self.num_sites)
            #return NotImplemented
        ham = copy(self.ham)
        ham += other.ham
        h = TIHamiltonian(self.num_sites, ham)
        h.cleanHam()
        return h

    def __sub__(self, other) :
        if self.num_sites != other.num_sites :
            self.extend_sites(other.num_sites)
            other.extend_sites(self.num_sites)
            #return NotImplemented
        ham = copy(self.ham)
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
        ham = copy(self.ham)
        for i in range(len(ham)) :
            ham[i] = (ham[i][0], ham[i][1] * other)
        h = TIHamiltonian(self.num_sites, ham)
        return h

    def __mul__(self, other) :
        # need more typing restrictions
        if type(other) == int  or  type(other) == float  or  type(other) == complex  or  isinstance(other, Expression) :
            return self.scalar_mul(other)
        if self.num_sites != other.num_sites :
            self.extend_sites(other.num_sites)
            other.extend_sites(self.num_sites)
            #return NotImplemented
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
        ham = copy(self.ham)
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

Empty = TIHamiltonian.empty()
