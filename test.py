from simuq.solver import generate_as

from systems.heis3 import qs
#from systems.single_x import qs

from aais.ibm7 import mach

generate_as(qs, mach, 5)
#generate_as(FloquetQS, FluxSCMach)
