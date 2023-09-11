from setuptools import find_namespace_packages, setup

with open("src/simuq/_version.py") as f:
    version = f.readlines()[-1].split()[-1].strip("\"'")

setup(
    name="SimuQ",
    version=version,
    description="A Python package for quantum simulation with analog compilation",
    url="https://github.com/PicksPeng/SimuQ",
    license="BSD 3-Clause Clear",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
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
        "test": ["tox"],
    },
)
