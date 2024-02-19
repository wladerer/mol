from pyscf import gto

def geometry_from_file(file: str, charge: int | None = None, spin: int | None = None, basis: int | None= None) -> gto.Mole:
    """
    Read the geometry from a file and return a gto.Mole object.

    Args:
        file (str): The path to the file containing the geometry.

    Returns:
        gto.Mole: The gto.Mole object representing the geometry.
    """
    mol = gto.Mole()
    mol.build(atom=file, charge=charge, spin=spin, basis=basis) # read the geometry from a file
    return mol