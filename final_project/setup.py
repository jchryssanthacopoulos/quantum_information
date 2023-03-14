"""Set up tebd package."""

from setuptools import find_packages
from setuptools import setup


setup(
    name="tebd",
    description="Implementation of the time evolving block decimation algorithm for matrix product states",
    url="https://github.com/jchryssanthacopoulos/quantum_information/tree/main/final_project",
    packages=find_packages(),
    install_requires=[
        "autoray==0.6.1",
        "jupyter==1.0.0",
        "numpy==1.23.5",
        "opt-einsum==3.3.0",
        "quimb==1.4.2",
        "scipy==1.10.1"
    ]
)
