import numpy as np
from braket.ahs import AnalogHamiltonianSimulation, AtomArrangement, DrivingField
from braket.aws import AwsDevice
from braket.timings.time_series import TimeSeries

from simuq import _version


def gen_braket_code(pos, clocks, pulse):
    # print(pos)
    # print(clocks)
    # # pulse = np.array(pulse)
    # print(pulse)
    # print(pulse.shape)
    import matplotlib.pyplot as plt

    register = AtomArrangement()
    for posi in pos:
        register.add(np.array([posi[0] * 1e-6, posi[1] * 1e-6]))

    delta = TimeSeries()
    omega = TimeSeries()
    phi = TimeSeries()

    delta.put(0, pulse[0][0])
    delta.put(clocks[1], pulse[0][0])
    for i in range(len(clocks) - 3):
        delta.put(clocks[i + 2], pulse[0][i])
    delta.put(clocks[-1], pulse[0][-1])
    omega.put(0, 0)
    omega.put(clocks[1], pulse[1][0])
    for i in range(len(clocks) - 3):
        omega.put(clocks[i + 2], pulse[1][i])
    omega.put(clocks[-1], 0)
    phi.put(0, 0)
    phi.put(clocks[1], pulse[2][0])
    for i in range(len(clocks) - 3):
        phi.put(clocks[i + 2], pulse[2][i])
    phi.put(clocks[-1], 0)

    drive = DrivingField(amplitude=omega, phase=phi, detuning=delta)

    ahs_program = AnalogHamiltonianSimulation(register=register, hamiltonian=drive)

    aquila_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/Aquila")
    user_agent = f"SimuQ/{_version.__version__}"
    aquila_qpu.aws_session.add_braket_user_agent(user_agent)
    discretized_ahs_program = ahs_program.discretize(aquila_qpu)

    return discretized_ahs_program


def gen_clocks(times, ramp_time=0.05):
    clocks = [0, ramp_time]
    for t in times:
        clocks.append(clocks[-1] + t)
    clocks.append(clocks[-1] + ramp_time)
    clocks = [t * 1e-6 for t in clocks]
    return clocks


def clean_as(alignment, sol_gvars, boxes):
    pos = [(0.0, 0.0)]
    for i in range(int(len(sol_gvars) / 2)):
        pos.append((sol_gvars[2 * i], sol_gvars[2 * i + 1]))
    # Fix rotation angle
    theta = -np.arctan2(pos[1][1], pos[1][0])
    print("Fixing rotation by ", theta)
    for i in range(len(pos)):
        new_x = np.cos(theta) * pos[i][0] - np.sin(theta) * pos[i][1]
        new_y = np.sin(theta) * pos[i][0] + np.cos(theta) * pos[i][1]
        pos[i] = (new_x, new_y)
    m = len(boxes)
    times = []
    pulse = [[0 for i in range(m)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for (i, j), ins, h, ins_lvars in box[0]:
            print(f"i={i}, ins_lvars={ins_lvars}")
            if i == 0:
                pulse[0][evo_idx] = ins_lvars[0] * 1e6
            else:
                pulse[1][evo_idx] = ins_lvars[0] * 1e6
                pulse[2][evo_idx] = ins_lvars[1] * 1e6
    return pos, gen_clocks(times), pulse


def transpile(alignment, sol_gvars, boxes, edges):
    # print(alignment)
    # print(sol_gvars)
    print(boxes)
    code = gen_braket_code(*clean_as(alignment, sol_gvars, boxes))
    return code
