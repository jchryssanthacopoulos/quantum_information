"""Set up tebd package."""

from setuptools import find_packages
from setuptools import setup


setup(
    name="tebd",
    description="Implementation of the time evolving block decimation algorithm for matrix product states",
    url="https://github.com/jchryssanthacopoulos/quantum_information/tree/main/final_project",
    packages=find_packages(),
    install_requires=[
        "numpy==1.23.5",
        "quimb==1.4.2",
        "scipy==1.10.1"
    ]
)
