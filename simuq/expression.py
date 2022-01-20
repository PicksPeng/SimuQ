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

    @classmethod
    def sin(cls, e) :
        if type(e) == int  or  type(e) == float  or  type(e) == complex :
            return math.sin(e)
        if isinstance(e, BaseVar) :
            e = e.to_exp()
        exp = lambda gvars, lvars : math.sin(e.exp(gvars, lvars))
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
