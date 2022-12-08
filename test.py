import time
from simuq.solver import generate_as

N = 20

qsys = 'MIS'
if qsys == 'MIS' :
    from systems.mis import GenQS
    qs = GenQS('Chain', k = N, dis_num = 1)
elif qsys == 'Heis' :
    from system.heis5 import qs

# from systems.heis5 import qs
# from systems.single_x import qs
# from systems.annealing import qs


#mode = 'IBM OpenPulse'
mode = 'Rigetti'
if mode == 'Circuit' :
    # Circuit / Ion Trap Setting
    from aais.iontrap import GenMach
    from backends.qiskit_iontrap import transpile
    n_mach = N
    start_time = time.time()
    mach = GenMach(n_mach)
    circ = transpile(n_mach, *generate_as(qs, mach, 1, 0.5))
    end_time = time.time()
    print(f'Time consumption: {end_time - start_time} s')
elif mode == 'Rydberg1d Bloqade' :
    # Bloqade Rydberg Setting
    from aais.rydberg1d import GenMach
    from backends.bloqade_rydberg import transpile
    n_mach = N
    start_time = time.time()
    mach = GenMach(n_mach)
    code = transpile(*generate_as(qs, mach, 1, 0.5))
    end_time = time.time()
    print(f'Time consumption: {end_time - start_time} s')
elif mode == 'Rydberg2d Bloqade' :
    # Bloqade Rydberg Setting
    from aais.rydberg2d import GenMach
    from backends.bloqade_rydberg import transpile
    n_mach = N
    start_time = time.time()
    mach = GenMach(n_mach)
    code = transpile(*generate_as(qs, mach, 1, 0.5))
    end_time = time.time()
    print(f'Time consumption: {end_time - start_time} s')
elif mode == 'Rydberg QuEra' :
    # Braket QuEra Rydberg Setting
    from aais.rydberg1d_aws import GenMach
    from backends.bloqade_rydberg_aws import transpile
    n_mach = N
    mach = GenMach(n_mach)
    code = transpile(*generate_as(qs, mach, 1))
elif mode == 'IBM OpenPulse' :
    # IBM Qiskit OpenPulse Setting
    from qiskit.providers.fake_provider import FakeGuadalupe, FakeToronto, FakeWashington

    from aais.ibm import get_mach
    from qiskit import IBMQ
    
    start_time = time.time()
    
    #backend=FakeGuadalupe()
    backend=FakeWashington()
    mach=get_mach(backend)
    assembly=generate_as(qs, mach, 1)
    with open("qaoa.as","w+") as f:
        f.write(str(assembly[1])+"\n")
        for item in assembly[2]:
            f.write(str(item)+"\n")
        for item in assembly[3]:
            f.write(str(item)+"\n")

    from Analog_Hamiltonian_Simulation.Analog_Hamiltonian_Simulator.IBM_Machine_new import IBM_Machine
    from Analog_Hamiltonian_Simulation.Analog_Hamiltonian_Simulator.Program import Program

    machine = IBM_Machine(backend)
    program = Program(machine)
    program.init_from_file("qaoa.as")

    #Schedule the instructions greedily using the sorted DAG
    program.schedule()

    # program.concrete_schedule.draw()
    # program.concrete_schedule.simulate()

    # generate the pulse schedule
    program.transpile()
    # program.pulse_schedule.draw()

    program.pulse_schedule.generate_external_schedule()
    # program.pulse_schedule.draw_external_schedule()
    
    end_time = time.time()
    print(f'Time consumption: {end_time - start_time} s')
elif mode == 'Rigetti' :
    from aais.rigetti import get_mach
    mach = get_mach()
    assembly = generate_as(qs, mach, 1)
    with open("qaoa_rigetti.as","w+") as f:
        f.write(str(assembly[1])+"\n")
        for item in assembly[2]:
            f.write(str(item)+"\n")
        for item in assembly[3]:
            f.write(str(item)+"\n")
    
    
# # IBM Qiskit 5qubit Setting
# '''
# from aais.ibm5 import mach
# from backends.qiskit_5 import transpile
# circ=transpile(*generate_as(qs, mach, 1, 0.5))
# '''



# # generate_as(FloquetQS, FluxSCMach)
# # generate_as(qs, mach, 1)


# # from aais.rydberg2d import Rydberg as mach
# # from backends.bloqade_rydberg2d import transpile


# #assembly=generate_as(qs, mach, 1)

# """
# circ=transpile(*generate_as(qs, mach, 1, 0.5))
# import pickle

# with open('circ.pickle', 'wb') as f:
#     pickle.dump(circ, f)
# """
