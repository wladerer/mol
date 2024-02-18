#!/usr/bin/env python3
from pyscf import gto
from pyscf import scf
from pyscf.geomopt.geometric_solver import optimize
# from pyscf import solvent

method_dict = {
    'rhf': scf.RHF,
    'uhf': scf.UHF,
    'rohf': scf.ROHF,
    'ghf': scf.GHF,
    'rks': scf.RKS,
    'uks': scf.UKS,
}

def hartee_to_eV(value: float) -> float:
    return value * 27.21138602

def bohr_to_angstrom(value: float) -> float:
    return value * 0.52917721067

def convert_units(args):
    if args.total_energy is not None:
        args.total_energy = hartee_to_eV(args.total_energy)
    if args.grms is not None:
        args.grms = hartee_to_eV(args.grms)
    if args.gmax is not None:
        args.gmax = hartee_to_eV(args.gmax)
    if args.drms is not None:
        args.drms = bohr_to_angstrom(args.drms)
    if args.dmax is not None:
        args.dmax = bohr_to_angstrom(args.dmax)

def geometry_from_file(file: str) -> gto.Mole:
    """
    Read the geometry from a file and return a gto.Mole object.

    Args:
        file (str): The path to the file containing the geometry.

    Returns:
        gto.Mole: The gto.Mole object representing the geometry.
    """
    mol = gto.Mole()
    mol.build(atom=file) # read the geometry from a file
    return mol

def geometry_opt(mol: gto.Mole, method: str, stdout=None, output=None) -> gto.Mole:
    """
    Perform geometry optimization on a molecular system.

    Args:
        mol (gto.Mole): The molecular system to optimize.
        method (str): The optimization method to use.
        stdout (file-like object, optional): The standard output stream for the optimization process.
        output (file-like object, optional): The output file to save the optimization results.

    Returns:
        gto.Mole: The optimized molecular system.
    """

    if stdout is not None:
        mol.stdout = stdout
    
    opt_mol = optimize(method_dict[method](mol), maxsteps=100)

    if output is not None:
        opt_mol.tofile(output, format='xyz')

    return opt_mol


def run(args):
    
    convert_units(args)

    mol = geometry_from_file(args.input)

    #update mol with user input
    if args.charge is not None:
        mol.charge = args.charge
    if args.spin is not None:
        mol.spin = args.spin
    if args.multiplicity is not None:
        mol.multiplicity = args.multiplicity
    if args.xc is not None:
        mol.xc = args.xc
    if args.basis is not None:
        mol.basis = args.basis
    
    
    geometry_opt(mol, args.method, args.stdout, args.output)

    return None