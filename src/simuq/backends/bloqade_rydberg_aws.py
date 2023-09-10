# def gen_bloqade_code(pos, clocks, pulse) :
#     n = len(pos)
#     code = "using Bloqade\n"
#     pos_1d = list(map(lambda x : (x,), pos))
#     code += f'atoms = AtomList({ pos_1d }) \n'
#     code_delta = "delta = ["
#     code_omega = "omega = ["
#     code_phi = "phi = ["
#     for i in range(n) :
#         code_delta += f'piecewise_constant(clocks = { clocks }, values = { pulse[0][i] }), '
#         code_omega += f'piecewise_constant(clocks = { clocks }, values = { pulse[1][i] }), '
#         code_phi += f'piecewise_constant(clocks = { clocks }, values = { pulse[2][i] }), '
#     code_omega += "]\n"
#     code_phi += "]\n"
#     code_delta += "]\n"
#     code += code_omega + code_phi + code_delta
#     code += "h = rydberg_h(atoms; Δ = delta, Ω = omega, ϕ = phi)"
#     return code

import numpy as np


def gen_braket_code(pos, clocks, pulse):
    # print(pos)
    # print(clocks)
    # # pulse = np.array(pulse)
    # print(pulse)
    # print(pulse.shape)
    from braket.ahs.atom_arrangement import AtomArrangement

    register = AtomArrangement()
    for posi in pos:
        register.add(np.array([0, posi]))
    from braket.ahs.driving_field import DrivingField
    from braket.ahs.time_series import TimeSeries

    delta = TimeSeries()
    omega = TimeSeries()
    phi = TimeSeries()
    # TODO piecewise constant or linear?

    delta.put(0, pulse[0][0][0])
    for i in range(len(clocks) - 1):
        delta.put((clocks[i] + clocks[i + 1]) / 2, pulse[0][0][i])
    delta.put(clocks[-1] / 2, pulse[0][0][-1])
    omega.put(0, pulse[0][1][0])
    for i in range(len(clocks) - 1):
        omega.put((clocks[i] + clocks[i + 1]) / 2, pulse[1][0][i])
    omega.put(clocks[-1] / 2, pulse[0][1][-1])
    phi.put(0, pulse[0][2][0])
    for i in range(len(clocks) - 1):
        phi.put((clocks[i] + clocks[i + 1]) / 2, pulse[2][0][i])
    phi.put(clocks[-1] / 2, pulse[0][2][-1])

    drive = DrivingField(amplitude=omega, phase=phi, detuning=delta)
    from braket.ahs.analog_hamiltonian_simulation import AnalogHamiltonianSimulation

    ahs_program = AnalogHamiltonianSimulation(register=register, hamiltonian=drive)
    from braket.aws import AwsDevice

    aquila_qpu = AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/Aquila")
    discretized_ahs_program = ahs_program.discretize(aquila_qpu)
    from braket.devices import LocalSimulator

    device = LocalSimulator("braket_ahs")

    result = device.run(discretized_ahs_program, shots=1_000_000).result()


def gen_clocks(times):
    clocks = [0]
    for t in times:
        clocks.append(clocks[-1] + t)
    return clocks


def clean_as(alignment, sol_gvars, boxes):
    pos = [0, *sol_gvars]
    n = len(pos)
    m = len(boxes)
    times = []
    pulse = [[0 for i in range(m)] for k in range(3)]
    for evo_idx in range(m):
        box = boxes[evo_idx]
        times.append(box[1])
        for (i, j), ins, h, ins_lvars in box[0]:
            # print(f"i={i}")
            if i == 0:
                pulse[0][evo_idx] = ins_lvars[0]
            else:
                pulse[1][evo_idx] = ins_lvars[0]
                pulse[2][evo_idx] = ins_lvars[1]

    return (pos, gen_clocks(times), pulse)


def transpile(alignment, sol_gvars, boxes, edges):
    code = gen_braket_code(*clean_as(alignment, sol_gvars, boxes))
    with open("transpiled.jl", "w") as f:
        print(code, file=f)
    return code


pulse = [
    [
        [
            -2.968997921691882,
            -2.302331255350557,
            -1.6356645886093095,
            -0.968997921317156,
            -0.30233125458228,
            0.36433541203769715,
            1.0310020784655602,
            1.6976687450141825,
            2.3643354106343843,
        ],
        [
            -3.000488222942375,
            -2.3338215563147067,
            -1.66715488935505,
            -1.0004882231196242,
            -0.33382155630970944,
            0.33284511042370446,
            0.9995117774412948,
            1.6661784439446585,
            2.3328451106648456,
        ],
        [
            -2.9689979249397465,
            -2.3023312582224844,
            -1.635664591334513,
            -0.9689979256061668,
            -0.30233125917951675,
            0.3643354077036973,
            1.0310020757343243,
            1.6976687415234875,
            2.3643354085635644,
        ],
        [
            -1.4708215063341783e-09,
            -1.428046216520107e-09,
            -1.4831853377199497e-09,
            -1.570594085422478e-09,
            -1.6746726602059686e-09,
            -2.5665215817997365e-09,
            -1.9309168757674367e-09,
            -1.8174553832496394e-09,
            -4.579114160373652e-09,
        ],
        [
            -1.8602264388175136e-09,
            -4.326363158835265e-10,
            0,
            -4.667673470364462e-10,
            -4.113543169861989e-10,
            -2.8613355627097614e-09,
            -4.031977288615272e-09,
            -1.1240022769670263e-09,
            -4.996166859307392e-10,
        ],
    ]
]
