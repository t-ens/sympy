from sympy.core import Basic, Expr, Mul, Add
from sympy.core.cache import cacheit
from sympy.utilities.iterables import variations 
from sympy.printing.str import *
from functools import lru_cache
from random import randint

class LyndonWord(Expr):
    """
    Provides basic operations on Lyndon words useful for Lie algebras.
    """

    is_commutative = False
    is_number = False

    def __new__(cls, word = ""):
        if not isinstance(word, str):
            breakpoint()
            #raise ValueError("Lyndon words must be strings")
        else:
            right_factors = [ word[i:] for i in range(1, len(word)) ]
            #A Lyndon word is lexicographically smaller than all its right factors
            if all([word < factor for factor in right_factors]):
                obj = Expr.__new__(cls, word)
                obj._word = word
                return obj
            else:
                raise ValueError("not a Lyndon word")

    @property
    def word(self):
        return self._word

    @property
    def deg(self):
        return len(self._word)

    def lyndon_factorization(self):
        """
        The Lyndon factorization of a Lyndon word w is the pair (w1,w2)
        where w2 is the minimal right factor of w and w1 is obtained by
        removing w2 from w
        """
        if self.deg <= 1: 
            return self
        else:
            right_factor = min([self.word[i:] for i in range(1, self.deg)])
            return LyndonWord(self.word[0:self.deg-len(right_factor)]), LyndonWord(right_factor)

    def bracket_form(self):
        """
        A Lyndon word corresponds to an element of the free Lie algebra by
        iteratively applying the Lie bracket to the pair obtained from the
        Lyndon factorization of a word.  
        """
        if self.deg <= 1:
            return self.word
        else:
            w1, w2 = self.lyndon_factorization()
            return "[" + w1.bracket_form() + "," + w2.bracket_form() + "]"

    def _sympystr(self, printer):
        return self.bracket_form()

    @property
    def lw_adjoint(self):
        """
        Given a Lie word w w.lw_adjoint is the function u = ad_w defined by u(z) = [w,z].
        """
        def u(z):
            if self == z:
                return 0
            elif z.word < self.word:
                return -bracket(z,self)
            elif self.deg == 1:
                return LyndonWord(self.word + z.word)
            else:
                x, y = self.lyndon_factorization()
                if y.word >= z.word:
                    return LyndonWord(self.word + z.word) 
                else:
                    return bracket(x, y.lw_adjoint(z)) + bracket(x.lw_adjoint(z), y)
        return u

    @staticmethod
    def generate(n, alph):
        """ 
        Returns a generator for  all Lyndon words of length at most n in the
        alphabet alph using Duval's algorithm. alph must be a string containing the
        letters of the desired alphabet in order from smallest to largest.
        """
        index=-1
        word = [alph[0]]
        while word:
            if index > -1:
                index = alph.index(word[-1]) + 1
            else:
                index = 0
            word[-1] = alph[index]
            yield "".join(word)
            m = len(word)
            while len(word) < n:
                word.append(word[-m])
            while word and word[-1] == alph[-1]:
                word.pop()

def bracket(v1, v2): 
    if isinstance(v1, Mul):
        return v1.args[0] * bracket(v1.args[1], v2)
    if isinstance(v2, Mul):
        return v2.args[0] * bracket(v1, v2.args[1])
    if isinstance(v1, Add):
        return sum([bracket(arg,v2) for arg in v1.args])
    if isinstance(v2, Add):
        return sum([bracket(v1,arg) for arg in v2.args])
    if isinstance(v1, LyndonWord) and isinstance(v2,LyndonWord):
        return v1.lw_adjoint(v2)

    raise ValueError("bracket only defined for linear combinations of Lyndon words")

class LieSeries(Expr):
    """
    Implement Lie series. A Lie series is defined by a function ser() which
    """

    show_degree = 3
    
    def __new__(cls, rule):
        obj = Expr.__new__(cls)
        obj._rule = rule
        obj._computed_vals = []
        return obj

    def ser(self, d):
      n = len(self._computed_vals)
      if len(self._computed_vals) < d:
          for x in range(n,d):  
              new_term = self._rule(x+1)
              if not isinstance(new_term, int):
                  self._computed_vals.append(new_term.expand())
              else:
                  self._computed_vals.append(new_term)
      return self._computed_vals[d-1]
      

    def __add__(self, other):
        return self.add(other) 

    def __mul__(self, c):
        return self.mul(c)

    def __rmul__(self, c):
        return self.mul(c)

    def add(self, other):
        return LieSeries(lambda d : self._rule(d) + other._rule(d))

    def mul(self, c):
        return LieSeries(lambda d : c*self._rule(d))

    def bracket(self, other):
        return LieSeries(lambda d : sum([bracket(self.ser(k), other.ser(d-k)) for k in range(1, d)]))

    def _sympystr(self, printer):
        #string = "".join([ printer._print(self.ser(i)) + " + " for i in range(1, self.show_degree+1) if not self.ser(i) == LyndonWord("")])
        #return string[:-3]
        return printer._print([self.ser(i) for i in range(1, self.show_degree+1)])

    @staticmethod
    def random(alph):
        def random_word_sum(d):
            words = [x for x  in LyndonWord.generate(d, alph) if len(x) == d]
            return sum([LyndonWord(words[randint(0, len(words)-1)]) for i in range(2,6)])
        return LieSeries(random_word_sum)

    @staticmethod
    def unknown(alph, name):
        def func(n):
            if n == "setter":
                return Null
            else:
                return
