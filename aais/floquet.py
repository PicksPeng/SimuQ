from simuq.environment import qubit, fock
from simuq.qmachine import *
from simuq.expression import Expression
from simuq.hamiltonian import TIHamiltonian
import numpy as np

# The Floquet Hamiltonian system in arXiv:2107.07311

L = 5
Js = [2*np.pi*0.01084, 2*np.pi*0.00028] # We assume they are transverse here.
l = len(Js)
geps = 0.1
hs = [0.1, 0.2, 0.3, 0.4, 0.5]

mach = QMachine()
qs = [qubit(mach) for i in range(L)]

Line = SignalLine(mach)

ins1 = Instruction(Line, 'native', 'H_Flip')
hflip = TIHamiltonian.empty(L)
for i in range(L) :
    hflip += geps / 2 * qs[i].X()
ins1.set_ham(hflip)

ins2 = Instruction(Line, 'native', 'H_Disorder')
hdis = TIHamiltonian.empty(L)
for i in range(L) :
    hdis += hs[i] / 2 * qs[i].Z()
ins2.set_ham(hdis)

ins3 = Instruction(Line, 'native', 'H_Int')
hint = TIHamiltonian.empty(L)
for j in range(1, l + 1) :
    for i in range(L - j) :
        hint += Js[j-1] / 2 * (qs[i].X() * qs[i+j].X() + qs[i].Y() * qs[i+j].Y())
ins3.set_ham(hint)
