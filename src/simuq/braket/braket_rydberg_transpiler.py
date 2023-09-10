import numpy as np
from braket.ahs import AnalogHamiltonianSimulation, AtomArrangement, DrivingField
from braket.timings.time_series import TimeSeries

from simuq.transpiler import Transpiler


def gen_braket_code(pos, clocks, pulse):
    # print(pos)
    # print(clocks)
    # # pulse = np.array(pulse)
    # print(pulse)
    # print(pulse.shape)

    register = AtomArrangement()
    for posi in pos:
        register.add(np.array([posi[0] * 1e-6, posi[1] * 1e-6]))

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

    ahs_program = AnalogHamiltonianSimulation(register=register, hamiltonian=drive)
    return ahs_program


def gen_clocks(times, ramp_time, state_prep_time):
    clocks = [0, ramp_time]
    for t in state_prep_time:
        clocks.append(clocks[-1] + t)
    clocks.append(clocks[-1] + ramp_time)
    for t in times:
        clocks.append(clocks[-1] + t)
    clocks.append(clocks[-1] + ramp_time)
    clocks = [t * 1e-6 for t in clocks]
    return clocks


def _initialize_positions_1d(sol_gvars):
    pos = [(0.0, 0.0)]
    for i in range(len(sol_gvars)):
        pos.append((sol_gvars[i], 0.0))
    return pos


def _initialize_positions_2d(sol_gvars, verbose):
    pos = [(0.0, 0.0)]
    for i in range(int(len(sol_gvars) / 2)):
        pos.append((sol_gvars[2 * i], sol_gvars[2 * i + 1]))
    # Fix rotation angle
    theta = -np.arctan2(pos[1][1], pos[1][0])
    if verbose > 0:
        print("Fixing rotation by ", theta)
    for i in range(len(pos)):
        new_x = np.cos(theta) * pos[i][0] - np.sin(theta) * pos[i][1]
        new_y = np.sin(theta) * pos[i][0] + np.cos(theta) * pos[i][1]
        pos[i] = (new_x, new_y)
    return pos


def clean_as(sol_gvars, boxes, ramp_time, state_prep, dimension, verbose=0):
    if dimension == 1:
        pos = _initialize_positions_1d(sol_gvars)
    elif dimension == 2:
        pos = _initialize_positions_2d(sol_gvars, verbose)
    else:
        raise NotImplementedError

    m = len(boxes)
    times = []
    pulse = [[0 for i in range(m)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for (i, j), ins, h, ins_lvars in box[0]:
            if verbose > 0:
                print(f"i={i}, ins_lvars={ins_lvars}")
            if i == 0:
                pulse[0][evo_idx] = ins_lvars[0] * 1e6
            else:
                pulse[1][evo_idx] = ins_lvars[0] * 1e6
                pulse[2][evo_idx] = ins_lvars[1] * 1e6
    if len(state_prep["times"]) > 0:
        pulse[0] = [v * 1e6 for v in state_prep["delta"]] + [pulse[0][0]] + pulse[0]
        pulse[1] = [v * 1e6 for v in state_prep["omega"]] + [pulse[1][0]] + pulse[1]
        pulse[2] = [v * 1e6 for v in state_prep["phi"]] + [pulse[2][0]] + pulse[2]
    else:
        for i in range(3):
            pulse[i] = [0, pulse[i][0]] + pulse[i]
    for i in range(len(pulse[0])):
        if pulse[1][i] < 0:
            pulse[1][i] = -pulse[1][i]
            pulse[2][i] += np.pi
    return pos, gen_clocks(times, ramp_time, state_prep["times"]), pulse


class BraketRydbergTranspiler(Transpiler):
    def __init__(self, dimension):
        self._dimension = dimension

    def transpile(self, sol_gvars, boxes, edges, ramp_time=0.05, state_prep=None, verbose=0):
        # print(sol_gvars)
        # print(boxes)
        code = gen_braket_code(
            *clean_as(sol_gvars, boxes, ramp_time, state_prep or {}, self._dimension, verbose)
        )
        return code
