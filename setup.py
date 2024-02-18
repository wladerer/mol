from setuptools import find_packages, setup

setup(
    name="mol",
    version="0.01",
    packages=find_packages(),
    install_requires=[
        "ase",
        "numpy",
        "matplotlib",
        "pymatgen",
        "mp-api",
        "pydantic",
    ],
    entry_points={"console_scripts": ["mol = mol.cli:main"]},
)
