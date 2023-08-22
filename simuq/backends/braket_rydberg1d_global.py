import numpy as np

def gen_braket_code(pos, clocks, pulse):
    # print(pos)
    # print(clocks)
    # # pulse = np.array(pulse)
    #print(pulse)
    # print(pulse.shape)
    import matplotlib.pyplot as plt
    from braket.ahs.atom_arrangement import AtomArrangement

    register = AtomArrangement()
    for posi in pos:
        register.add(np.array([posi[0]*1e-6, posi[1]*1e-6]))

    from braket.ahs.driving_field import DrivingField
    from braket.timings.time_series import TimeSeries

    delta = TimeSeries()
    omega = TimeSeries()
    phi = TimeSeries()

    delta.put(0, pulse[0][0])
    for i in range(len(clocks) - 2):
        delta.put(clocks[i + 1], pulse[0][i])
    delta.put(clocks[-1], pulse[0][-1])
    omega.put(0, 0)
    for i in range(len(clocks) - 2):
        omega.put(clocks[i + 1], pulse[1][i])
    omega.put(clocks[-1], 0)
    phi.put(0, 0)
    for i in range(len(clocks) - 2):
        phi.put(clocks[i + 1], pulse[2][i])
    phi.put(clocks[-1], 0)

    drive = DrivingField(amplitude=omega, phase=phi, detuning=delta)
    from braket.ahs.analog_hamiltonian_simulation import AnalogHamiltonianSimulation

    ahs_program = AnalogHamiltonianSimulation(register=register, hamiltonian=drive)

    try :
        from braket.aws import AwsDevice
        aquila_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/Aquila")
    except :
        raise Exception("Please check your configurations of IAM roles, regions, and credentials according to instructions in: https://boto3.amazonaws.com/v1/documentation/api/latest/guide/quickstart.html#configuration")

    discretized_ahs_program = ahs_program.discretize(aquila_qpu)

    return discretized_ahs_program
    
    from braket.devices import LocalSimulator
    device = LocalSimulator("braket_ahs")

    result = device.run(discretized_ahs_program, shots=1_000_000).result()


def gen_clocks(times, ramp_time, state_prep_time):
    clocks = [0, ramp_time]
    for t in state_prep_time :
        clocks.append(clocks[-1] + t)
    clocks.append(clocks[-1] + ramp_time)
    for t in times:
        clocks.append(clocks[-1] + t)
    clocks.append(clocks[-1] + ramp_time)
    clocks = [t * 1e-6 for t in clocks]
    return clocks


def clean_as(sol_gvars, boxes, ramp_time, state_prep, verbose = 0):
    pos = [(0., 0.) ]
    for i in range(len(sol_gvars)) :
        pos.append((sol_gvars[i], 0.))
    n = len(pos)
    m = len(boxes)
    times = []
    pulse = [[0 for i in range(m)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for ((i, j), ins, h, ins_lvars) in box[0]:
            if verbose > 0 :
                print(f"i={i}, ins_lvars={ins_lvars}")
            if i == 0:
                pulse[0][evo_idx] = ins_lvars[0] * 1e6
            else:
                pulse[1][evo_idx] = ins_lvars[0] * 1e6
                pulse[2][evo_idx] = ins_lvars[1] * 1e6
    if len(state_prep['times']) > 0 :
        pulse[0] = [v * 1e6 for v in state_prep['delta']] + [pulse[0][0]] + pulse[0]
        pulse[1] = [v * 1e6 for v in state_prep['omega']] + [pulse[1][0]] + pulse[1]
        pulse[2] = [v * 1e6 for v in state_prep['phi']] + [pulse[2][0]] + pulse[2]
    else :
        for i in range(3) :
            pulse[i] = [0, pulse[i][0]] + pulse[i]
    for i in range(len(pulse[0])) :
        if pulse[1][i] < 0 :
            pulse[1][i] = -pulse[1][i]
            pulse[2][i] += np.pi
    return (pos, gen_clocks(times, ramp_time, state_prep['times']), pulse)


def transpile(sol_gvars, boxes, edges, ramp_time = 0.05, state_prep = {}, verbose = 0):
    #print(sol_gvars)
    #print(boxes)
    code = gen_braket_code(*clean_as(sol_gvars, boxes, ramp_time, state_prep, verbose))
    return code
