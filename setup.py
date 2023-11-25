from setuptools import find_namespace_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("src/simuq/_version.py") as f:
    version = f.readlines()[-1].split()[-1].strip("\"'")

setup(
    name="simuq",
    version=version,
    description="A Python package for quantum simulation with analog compilation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pickspeng.github.io/SimuQ/",
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
        "qiskit": [
            "matplotlib",
            "qiskit",
            "qiskit-terra",
            "qiskit-ibmq-provider",
        ],
        "qutip": [
            "qutip",
        ],
        "all": ["simuq[braket, ionq, qiskit, qutip, dev]"],
        "dev": ["tox"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
