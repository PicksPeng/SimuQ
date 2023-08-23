import codecs
import os.path

from setuptools import find_packages, setup


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="SimuQ",
    version=get_version("simuq/__init__.py"),
    description="A Python package for quantum simulation with analog compilation",
    url="https://github.com/PicksPeng/SimuQ",
    license="BSD 3-Clause Clear",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "networkx",
    ],
    extras_require={
        "braket": [
            "matplotlib",
            "amazon-braket-sdk",
        ],
        "ionq": [
            "requests",
        ],
        "dreal": [
            "dreal",
        ],
        "qiskit": [
            "matplotlib",
            "qiskit",
            "qiskit-terra",
        ],
        "qutip": [
            "qutip",
        ],
        "all": ["simuq[braket, ionq, dreal, qiskit, qutip]"],
        "dev": ["tox"],
    },
)
