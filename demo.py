import time
from simuq.solver import generate_as

K = 5
#N = K * K
N = K

Repetition = 1

D = 1

Trot = 1


qsys = 'MIS'

if qsys == 'MIS' :
    from systems.mis import GenQS
    qs = GenQS('Chain', k = K, dis_num = D)
elif qsys == 'Heis' :
    from systems.heis7 import qs
elif qsys == 'QAOA' :
    from systems.qaoa import GenQS
    qs = GenQS(n = K, p = D)



mode = 'IBM OpenPulse'

if mode == 'Circuit' :
    # Circuit / Ion Trap Setting
    from aais.iontrap import GenMach
    from backends.qiskit_iontrap import transpile
    n_mach = N
    start_time = time.time()
    for r in range(Repetition) :
        mach = GenMach(n_mach)
        circ = transpile(n_mach, *generate_as(qs, mach, Trot, 0.5))
        print(f'finishing {r}')
    end_time = time.time()
    print(f'Avg Time consumption: {(end_time - start_time) / Repetition} s')


elif mode == 'Rydberg1d Bloqade' :
    # Bloqade Rydberg Setting
    from aais.rydberg1d import GenMach
    from backends.bloqade_rydberg import transpile
    n_mach = N
    start_time = time.time()
    for r in range(Repetition) :
        mach = GenMach(n_mach)
        code = transpile(*generate_as(qs, mach, Trot, 0.5))
        print(f'finishing {r}')
    end_time = time.time()
    print(f'Avg Time consumption: {(end_time - start_time) / Repetition} s')


elif mode == 'Rydberg2d Bloqade' :
    # Bloqade Rydberg Setting
    from aais.rydberg2d import GenMach
    from backends.bloqade_rydberg2d import transpile
    n_mach = N
    start_time = time.time()
    for r in range(Repetition) :
        mach = GenMach(n_mach)
        code = transpile(*generate_as(qs, mach, Trot, 0.5))
        print(f'finishing {r}')
    end_time = time.time()
    print(f'Avg Time consumption: {(end_time - start_time) / Repetition} s')


elif mode == 'Rydberg QuEra' :
    # Braket QuEra Rydberg Setting
    from aais.rydberg1d_aws import GenMach
    from backends.bloqade_rydberg_aws import transpile
    n_mach = N
    mach = GenMach(n_mach)
    code = transpile(*generate_as(qs, mach, Trot))


elif mode == 'IBM OpenPulse' :
    # IBM Qiskit OpenPulse Setting
    from qiskit.providers.fake_provider import FakeJakarta, FakeGuadalupe, FakeToronto, FakeWashington
    backend=FakeJakarta()
    #backend=FakeGuadalupe()
    #backend=FakeToronto()
    #backend=FakeWashington()

    from aais.ibm import get_mach
    from backends.qiskit_pulse_ibm import transpile
    
    start_time = time.time()

    for r in range(Repetition) :
        mach = get_mach(backend)

        schedule = transpile(backend, *generate_as(qs, mach, Trot))
        
        print(f'finishing {r}')
    
    end_time = time.time()
    print(f'Avg Time consumption: {(end_time - start_time) / Repetition} s')

    ''' old version 
    # IBM Qiskit OpenPulse Setting
    from qiskit.providers.fake_provider import FakeJakarta, FakeGuadalupe, FakeToronto, FakeWashington

    from aais.ibm import get_mach
    from qiskit import IBMQ
    from backends.Analog_Hamiltonian_Simulation.Analog_Hamiltonian_Simulator.IBM_Machine_new import IBM_Machine
    from backends.Analog_Hamiltonian_Simulation.Analog_Hamiltonian_Simulator.Program import Program
    
    start_time = time.time()

    for r in range(Repetition) :

        #backend=FakeJakarta()
        backend=FakeGuadalupe()
        #backend=FakeToronto()
        #backend=FakeWashington()

        mach=get_mach(backend)

        assembly=generate_as(qs, mach, 1)
        
        program = Program(IBM_Machine(backend))
        program.init_from_as(assembly)

        #Schedule the instructions greedily using the sorted DAG
        program.schedule()

        # program.concrete_schedule.draw()
        # program.concrete_schedule.simulate()

        # generate the pulse schedule
        program.transpile()
        # program.pulse_schedule.draw()

        program.pulse_schedule.generate_external_schedule()
        # program.pulse_schedule.draw_external_schedule()
        
        print(f'finishing {r}')
    
    end_time = time.time()
    print(f'Avg Time consumption: {(end_time - start_time) / Repetition} s')
    '''


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
