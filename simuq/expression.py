'''
The class file describing expressions with variables.

Expressions are effectively functions taking a valuation of
variables and generating (real) results. We support the
natural construction of expressions, where users write 
+, -, *, / naturally on expressions.

These expressions will be used in representing the 
coefficients in machine's instruction Hamiltonians, where 
variables belonging to the instruction have effects on
the Hamiltonian.

The basic operations are overloaded.
'''

import math
import numpy as np

class BaseVar :
    """ The basic variables.

    The constrainsts specific to the variables are stored here,
    like the initial value, lower and upper bounds of them.
    """
    def __init__(self, mach) :
        self.mach = mach
        self.init_value = 0
        self.lower_bound = -np.inf
        self.upper_bound = np.inf

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
    """ The expressions.
    
    It is effectively a function taking a valuation of global variables
    and local variables and generating a number.
    """
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
        if type(e) == int  or  type(e) == float :
            return math.cos(e)
        if isinstance(e, BaseVar) :
            e = e.to_exp()
        if isinstance(e, np.float64) :
            return np.cos(e)

        def exp(gvars, lvars) :
            sub = e.exp(gvars, lvars)
            if type(sub) == int  or  type(sub) == float :
                return math.cos(sub)
            if isinstance(sub, np.float64) :
                return np.cos(sub)
            import dreal as dr
            if isinstance(sub, dr._dreal_py.Variable)  or  isinstance(sub, dr._dreal_py.Expression) :
                return dr.cos(sub)
            pp = qq
            return NotImplemented

        return cls(e.mach, exp)

    @classmethod
    def sin(cls, e) :
        if type(e) == int  or  type(e) == float :
            return math.sin(e)
        if isinstance(e, BaseVar) :
            e = e.to_exp()
        if isinstance(e, np.float64) :
            return np.sin(e)
        
        def exp(gvars, lvars) :
            sub = e.exp(gvars, lvars)
            if type(sub) == int  or  type(sub) == float :
                return math.cos(sub)
            if isinstance(sub, np.float64) :
                return np.sin(sub)
            import dreal as dr
            if type(sub) == dr._dreal_py.Variable  or  type(sub) == d._dreal_py.Expression :
                return dr.sin(sub)
            pp = qq
            return NotImplemented
        
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
            exp = lambda gvars, lvars : other / self.exp(gvars, lvars)
            e = Expression(self.mach, exp)
            return e
        if isinstance(other, BaseVar) :
            other = other.to_exp()
        if not hasattr(other, "exp") :
            return NotImplemented
        exp = lambda gvars, lvars : other.exp(gvars, lvars) / self.exp(gvars, lvars)
        e = Expression(self.mach, exp)
        return e

    def exp_eval(self, gvars, lvars) :
        return self.exp(gvars, lvars)
