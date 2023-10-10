from math import erf

import numpy as np
from qiskit import QuantumCircuit
from qiskit.pulse import (
    ControlChannel,
    DriveChannel,
    GaussianSquare,
    Play,
    Waveform,
    builder,
)
import enum



import numpy as np
from qiskit import QuantumCircuit, schedule, transpile
from qiskit.circuit import Gate
from math import erf
from typing import List, Tuple

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit.pulse import (
    ControlChannel,
    DriveChannel,
    GaussianSquare,
    Play,
    Waveform,
    builder,
)
from qiskit.pulse.filters import filter_instructions
from qiskit.pulse.instruction_schedule_map import InstructionScheduleMap


class CRCalType(enum.Enum):
    """Estimated calibration type of backend cross resonance operations."""

    ECR_FORWARD = "Echoed Cross Resonance corresponding to native operation"
    ECR_REVERSE = "Echoed Cross Resonance reverse of native operation"
    ECR_CX_FORWARD = "Echoed Cross Resonance CX corresponding to native operation"
    ECR_CX_REVERSE = "Echoed Cross Resonance CX reverse of native operation"
    DIRECT_CX_FORWARD = "Direct CX corresponding to native operation"
    DIRECT_CX_REVERSE = "Direct CX reverse of native operation"


def _check_calibration_type(
    inst_sched_map: InstructionScheduleMap, qubits: List[int]
) -> Tuple[CRCalType, List[Play], List[Play]]:
    """A helper function to check type of CR calibration.

    Args:
        inst_sched_map: instruction schedule map of the backends
        qubits: ordered tuple of qubits for cross resonance (q_control, q_target)

    Returns:
        Filtered instructions and most-likely type of calibration.

    Raises:
        QiskitError: Unknown calibration type is detected.
    """
    cal_type = None
    if inst_sched_map.has("cx", qubits):
        cr_sched = inst_sched_map.get("cx", qubits=qubits)
    elif inst_sched_map.has("ecr", qubits):
        cr_sched = inst_sched_map.get("ecr", qubits=qubits)
        cal_type = CRCalType.ECR_FORWARD
    elif inst_sched_map.has("ecr", tuple(reversed(qubits))):
        cr_sched = inst_sched_map.get("ecr", tuple(reversed(qubits)))
        cal_type = CRCalType.ECR_REVERSE
    else:
        raise QiskitError(
            f"Native direction cannot be determined: operation on qubits {qubits} "
            f"for the following instruction schedule map:\n{inst_sched_map}"
        )

    cr_tones = [t[1] for t in filter_instructions(cr_sched, [_filter_cr_tone]).instructions]
    comp_tones = [t[1] for t in filter_instructions(cr_sched, [_filter_comp_tone]).instructions]

    if cal_type is None:
        if len(comp_tones) == 0:
            raise QiskitError(
                f"{repr(cr_sched)} has no target compensation tones. "
                "Native ECR direction cannot be determined."
            )
        # Determine native direction, assuming only single drive channel per qubit.
        # This guarantees channel and qubit index equality.
        if comp_tones[0].channel.index == qubits[1]:
            cal_type = CRCalType.ECR_CX_FORWARD
        else:
            cal_type = CRCalType.ECR_CX_REVERSE

    if len(cr_tones) == 2 and len(comp_tones) in (0, 2):
        # ECR can be implemented without compensation tone at price of lower fidelity.
        # Remarkable noisy terms are usually eliminated by echo.
        return cal_type, cr_tones, comp_tones

    if len(cr_tones) == 1 and len(comp_tones) == 1:
        # Direct CX must have compensation tone on target qubit.
        # Otherwise, it cannot eliminate IX interaction.
        if comp_tones[0].channel.index == qubits[1]:
            return CRCalType.DIRECT_CX_FORWARD, cr_tones, comp_tones
        else:
            return CRCalType.DIRECT_CX_REVERSE, cr_tones, comp_tones
    raise QiskitError(
        f"{repr(cr_sched)} is undefined pulse sequence. "
        "Check if this is a calibration for cross resonance operation."
    )


def test_native(q, backend):
    inst_sched_map = backend.defaults().instruction_schedule_map
    a = _check_calibration_type(inst_sched_map, q)
    return a[0] in [CRCalType.ECR_CX_FORWARD, CRCalType.ECR_FORWARD]


def get_coupling_map_x(backend):
    coupling_map = backend.configuration().coupling_map
    n_qubits = backend.configuration().n_qubits
    coupling_dict = [[] for i in range(n_qubits)]
    coupling_dict_with_x = [[] for i in range(n_qubits)]

    for q in coupling_map:
        if q[0] not in coupling_dict[q[1]]:
            coupling_dict[q[1]].append(q[0])
        if q[1] not in coupling_dict[q[0]]:
            coupling_dict[q[0]].append(q[1])
    for i in range(n_qubits):
        for j in coupling_dict[i]:
            if test_native([j, i], backend):
                coupling_dict_with_x[i].append(j)
    return coupling_dict_with_x


def calc_area(duration, width, sigma):
    risefall_sigma_ratio = (duration - width) / sigma

    return width + sigma * np.sqrt(2 * np.pi) * erf(risefall_sigma_ratio)


def rescale_cr_inst(instruction: Play, theta: float, sample_mult: int = 16):
    params = instruction.pulse.parameters.copy()
    risefall_sigma_ratio = (params["duration"] - params["width"]) / params["sigma"]

    # The error function is used because the Gaussian may have chopped tails.
    # Area is normalized by amplitude.
    # This makes widths robust to the rounding error.
    risefall_area = params["sigma"] * np.sqrt(2 * np.pi) * erf(risefall_sigma_ratio)
    full_area = params["width"] + risefall_area

    # Get estimate of target area. Assume this is pi/2 controlled rotation.
    cal_angle = np.pi / 2
    target_area = abs(theta) / cal_angle * full_area
    new_width = target_area - risefall_area

    if new_width >= 0:
        width = new_width
        params["amp"] *= np.sign(theta)
    else:
        width = 0
        params["amp"] *= np.sign(theta) * target_area / risefall_area

    round_duration = (
        round((width + risefall_sigma_ratio * params["sigma"]) / sample_mult) * sample_mult
    )
    params["duration"] = round_duration
    params["width"] = width

    stretched_pulse = GaussianSquare(**params)
    return stretched_pulse, instruction.channel


def zxI_Ixz_pulse(qz_1, qz_2, qx, theta, backend):
    instruction_schedule_map = backend.defaults().instruction_schedule_map

    qubit_pair = [qz_1, qx]
    cal_type, cr_tones, comp_tones = _check_calibration_type(instruction_schedule_map, qubit_pair)
    zx1_cr_pulse, zx1_cr_chan = rescale_cr_inst(cr_tones[0], theta)
    zx1_comp_pulse, zx1_comp_chan = rescale_cr_inst(comp_tones[0], theta)

    qubit_pair = [qz_2, qx]
    cal_type, cr_tones, comp_tones = _check_calibration_type(instruction_schedule_map, qubit_pair)
    zx2_cr_pulse, zx2_cr_chan = rescale_cr_inst(cr_tones[0], theta)
    zx2_comp_pulse, zx2_comp_chan = rescale_cr_inst(comp_tones[0], theta)

    assert zx1_comp_chan == zx2_comp_chan
    if zx1_cr_pulse.duration < zx2_cr_pulse.duration:
        print(zx1_cr_pulse.duration)

        area_1 = calc_area(zx1_cr_pulse.duration, zx1_cr_pulse.width, zx1_cr_pulse.sigma)
        area_2 = calc_area(zx2_cr_pulse.duration, zx2_cr_pulse.width, zx2_cr_pulse.sigma)
        amp_factor = area_1 / area_2
        zx1_cr_pulse = GaussianSquare(
            duration=zx2_cr_pulse.duration,
            amp=zx1_cr_pulse.amp * amp_factor,
            sigma=zx1_cr_pulse.sigma,
            width=zx2_cr_pulse.width,
            angle=zx1_cr_pulse.angle,
        )
        zx1_comp_pulse = GaussianSquare(
            duration=zx2_cr_pulse.duration,
            amp=zx1_comp_pulse.amp * amp_factor,
            sigma=zx1_comp_pulse.sigma,
            width=zx2_comp_pulse.width,
            angle=zx1_comp_pulse.angle,
        )

    elif zx1_cr_pulse.duration > zx2_cr_pulse.duration:
        area_1 = calc_area(zx1_cr_pulse.duration, zx1_cr_pulse.width, zx1_cr_pulse.sigma)
        area_2 = calc_area(zx2_cr_pulse.duration, zx2_cr_pulse.width, zx2_cr_pulse.sigma)
        amp_factor = area_2 / area_1
        zx2_cr_pulse = GaussianSquare(
            duration=zx1_cr_pulse.duration,
            amp=zx2_cr_pulse.amp * amp_factor,
            sigma=zx2_cr_pulse.sigma,
            width=zx1_cr_pulse.width,
            angle=zx2_cr_pulse.angle,
        )
        zx2_comp_pulse = GaussianSquare(
            duration=zx1_cr_pulse.duration,
            amp=zx2_comp_pulse.amp * amp_factor,
            sigma=zx2_comp_pulse.sigma,
            width=zx1_comp_pulse.width,
            angle=zx2_comp_pulse.angle,
        )
    assert zx1_comp_pulse.sigma == zx2_comp_pulse.sigma
    assert zx1_comp_pulse.width == zx2_comp_pulse.width
    complex_amp=zx1_comp_pulse.amp*np.exp(1j*zx1_comp_pulse.angle)+zx2_comp_pulse.amp*np.exp(1j*zx2_comp_pulse.angle)
    new_amp=np.abs(complex_amp)
    new_angle=np.angle(complex_amp)
    zx12_comp_pulse = GaussianSquare(
        duration=zx1_cr_pulse.duration,
        amp=new_amp,
        sigma=zx1_comp_pulse.sigma,
        width=zx1_comp_pulse.width,
        angle=new_angle,
    )


    with builder.build() as sched:
        builder.play(zx1_cr_pulse, zx1_cr_chan)
        builder.play(zx2_cr_pulse, zx2_cr_chan)
        builder.play(zx12_comp_pulse, zx1_comp_chan)

    return sched


def _filter_cr_tone(time_inst_tup):
    """A helper function to filter pulses on control channels."""
    valid_types = ["GaussianSquare"]

    _, inst = time_inst_tup
    if isinstance(inst, Play) and isinstance(inst.channel, ControlChannel):
        pulse = inst.pulse
        if isinstance(pulse, Waveform) or pulse.pulse_type in valid_types:
            return True
    return False


def _filter_comp_tone(time_inst_tup):
    """A helper function to filter pulses on drive channels."""
    valid_types = ["GaussianSquare"]

    _, inst = time_inst_tup
    if isinstance(inst, Play) and isinstance(inst.channel, DriveChannel):
        pulse = inst.pulse
        if isinstance(pulse, Waveform) or pulse.pulse_type in valid_types:
            return True
    return False


def build_zxi_ixz(backend, with_measure = True):
    if isinstance(backend.name, str) :
        name = backend.name
    else :
        name = backend.name()
    if name == "ibmq_lima":
        qubit_z1 = 0
        qubit_x = 1
        qubit_z2 = 2
    elif name == "ibmq_guadalupe":
        qubit_z1 = 0
        qubit_x = 1
        qubit_z2 = 2
    elif name == "ibmq_jakarta":
        qubit_z1 = 0
        qubit_x = 1
        qubit_z2 = 2
    elif name == "ibm_lagos":
        qubit_z1 = 0
        qubit_x = 1
        qubit_z2 = 2
    elif name == "ibmq_mumbai":
        qubit_z1 = 0
        qubit_x = 1
        qubit_z2 = 2
    elif name == "ibm_cairo":
        qubit_z1 = 1
        qubit_x = 2
        qubit_z2 = 3
    elif name == "ibm_hanoi":
        qubit_z1 = 10
        qubit_x = 12
        qubit_z2 = 13
    else:
        raise ValueError("Unsupported backend")
    n_qubits = backend.configuration().n_qubits
    zxI_Ixz_sched = zxI_Ixz_pulse(qubit_z1, qubit_z2, qubit_x, 2, backend)
    zxI_Ixz_gate = Gate("zxI_Ixz_gate", 3, [])
    measure_gate = Gate("measure_gate", 3, [])
    circ = QuantumCircuit(n_qubits, 3)
    circ.append(zxI_Ixz_gate, [qubit_z1, qubit_z2, qubit_x])


    circ.add_calibration("zxI_Ixz_gate", [qubit_z1, qubit_z2, qubit_x], zxI_Ixz_sched)

    # circ.add_calibration("measure_gate", [qubit_z1, qubit_z2, qubit_x], measure_sched)
    circ = transpile(circ, backend)
    # circ.measure(qubit_z1,0)
    # circ.measure(qubit_z2,1)
    # circ.measure(qubit_x,2)
    if with_measure :
        circ.measure_all()
    # sched = schedule(circ, backend)
    # return sched.duration
    return circ