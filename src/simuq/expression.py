"""
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
"""

import math

import numpy as np


class BaseVar:
    """The basic variables.

    The constraints specific to the variables are stored here,
    like the initial value, lower and upper bounds of them.
    """

    def __init__(self, init_value=0, lower_bound=-np.inf, upper_bound=np.inf):
        self.init_value = init_value
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.index = None

    def to_exp(self):
        return Expression.id_var(self)

    # Index should allow the value of the BaseVar to be obtained from a list of values of global and local variables.
    # Global variables have indices from 0 to num_gvars - 1
    # Local variables have indices from num_gvars to num_gvars + num_lvars - 1
    def set_index(self, index):
        self.index = index

    def __neg__(self):
        return -self.to_exp()

    def __add__(self, other):
        return self.to_exp() + other

    def __mul__(self, other):
        return self.to_exp() * other

    def __sub__(self, other):
        return self.to_exp() - other

    def __truediv__(self, other):
        return self.to_exp() / other


# Given two lists, return the union, and the indices of the items of the original lists in the union.
def find_union_indices(list1, list2):
    union_list = list(set(list1).union(list2))
    union_dict = dict((item, index) for index, item in enumerate(union_list))
    indices1 = [union_dict[item] for item in list1]
    indices2 = [union_dict[item] for item in list2]
    return union_list, indices1, indices2


class Expression:
    """The expressions.

    It is effectively a function taking a valuation of global variables
    and local variables and generating a number.
    """

    def __init__(self, exp, vars):
        self.exp = exp
        self.vars = vars

    @classmethod
    def unit(cls):
        exp = lambda values: 1
        return cls(exp, [])

    @classmethod
    def id_var(cls, var):
        exp = lambda values: values[0]
        return cls(exp, [var])

    @classmethod
    def cos(cls, e):
        if isinstance(e, (int, float)):
            return math.cos(e)
        if isinstance(e, BaseVar):
            e = e.to_exp()
        if isinstance(e, np.float64):
            return np.cos(e)

        def exp(values):
            sub = e.exp(values)
            if isinstance(sub, (int, float)):
                return math.cos(sub)
            if isinstance(sub, np.float64):
                return np.cos(sub)
            return NotImplemented

        return cls(exp, e.vars)

    @classmethod
    def sin(cls, e):
        if isinstance(e, (int, float)):
            return math.sin(e)
        if isinstance(e, BaseVar):
            e = e.to_exp()
        if isinstance(e, np.float64):
            return np.sin(e)

        def exp(values):
            sub = e.exp(values)
            if isinstance(sub, (int, float)):
                return math.sin(sub)
            if isinstance(sub, np.float64):
                return np.sin(sub)
            return NotImplemented

        return cls(exp, e.vars)

    def __neg__(self):
        exp = lambda vars: -self.exp(vars)
        e = Expression(exp, self.vars)
        return e

    def __add__(self, other):
        if isinstance(other, (int, float, complex)):
            exp = lambda values: self.exp(values) + other
            e = Expression(exp, self.vars)
            return e
        if isinstance(other, BaseVar):
            other = other.to_exp()
        if not hasattr(other, "exp"):
            return NotImplemented

        new_vars, indices1, indices2 = find_union_indices(self.vars, other.vars)
        exp = lambda values: self.exp([values[index] for index in indices1]) + other.exp(
            [values[index] for index in indices2]
        )
        e = Expression(exp, new_vars)
        return e

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __radd__(self, other):
        if isinstance(other, (int, float, complex)):
            return self.__add__(other)
        elif isinstance(other, BaseVar):
            other = other.to_exp()
            return self.__add__(other)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float, complex)):
            exp = lambda values: self.exp(values) * other
            e = Expression(exp, self.vars)
            return e
        if isinstance(other, BaseVar):
            other = other.to_exp()
        if not hasattr(other, "exp"):
            return NotImplemented

        new_vars, indices1, indices2 = find_union_indices(self.vars, other.vars)
        exp = lambda values: self.exp([values[index] for index in indices1]) * other.exp(
            [values[index] for index in indices2]
        )
        e = Expression(exp, new_vars)
        return e

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex)):
            return self.__mul__(other)
        elif isinstance(other, BaseVar):
            other = other.to_exp()
            return self.__mul__(other)
        else:
            return NotImplemented

    def __pow__(self, other):
        if isinstance(other, (int, float, complex)):
            exp = lambda values: self.exp(values) ** other
            e = Expression(exp, self.vars)
            return e
        if isinstance(other, BaseVar):
            other = other.to_exp()
        if not hasattr(other, "exp"):
            return NotImplemented

        new_vars, indices1, indices2 = find_union_indices(self.vars, other.vars)
        exp = lambda values: self.exp([values[index] for index in indices1]) ** other.exp(
            [values[index] for index in indices2]
        )
        e = Expression(exp, new_vars)
        return e

    def __truediv__(self, other):
        if isinstance(other, (int, float, complex)):
            exp = lambda values: self.exp(values) / other
            e = Expression(exp, self.vars)
            return e
        if isinstance(other, BaseVar):
            other = other.to_exp()
        if not hasattr(other, "exp"):
            return NotImplemented
        new_vars, indices1, indices2 = find_union_indices(self.vars, other.vars)
        exp = lambda values: self.exp([values[index] for index in indices1]) / other.exp(
            [values[index] for index in indices2]
        )
        e = Expression(exp, new_vars)
        return e

    def __rtruediv__(self, other):
        if isinstance(other, (int, float, complex)):
            exp = lambda values: other / self.exp(values)
            e = Expression(exp, self.vars)
            return e
        if isinstance(other, BaseVar):
            other = other.to_exp()
        if not hasattr(other, "exp"):
            return NotImplemented

        new_vars, indices1, indices2 = find_union_indices(other.vars, self.vars)
        exp = lambda values: other.exp([values[index] for index in indices1]) / self.exp(
            [values[index] for index in indices2]
        )
        e = Expression(exp, new_vars)

        return e

    def exp_eval(self, gvars, lvars):
        # The weird comprehension below is to handle np array shapes
        values = [gvar for gvar in gvars] + [lvar for lvar in lvars]
        return_val = self.exp([values[var.index] for var in self.vars])
        return return_val
