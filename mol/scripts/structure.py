#!/usr/bin/env python3

import numpy as np
from ase.io import read
import json


def get_atoms(args):
    """Creates ASE atoms object from a file"""

    atoms = read(args.input)

    return atoms


def sort_poscar(args):
    """Sorts a POSCAR file"""
    from pymatgen.core import Structure
    from pymatgen.io.vasp.inputs import Poscar

    structure = Structure.from_file(args.input)
    poscar = Poscar(structure, sort_structure=True)

    if not args.output:
        print(poscar.get_str())
    else:
        poscar.write_file(f"{args.output}")

    return poscar


def structure_from_mpi_code(mpcode: str, api_key: str, is_conventional: bool = True):
    """
    Creates a pymatgen structure from a code
    """
    from mp_api.client import MPRester

    if not mpcode.startswith("mp-"):
        mpcode = "mp-" + mpcode

    with MPRester(api_key, mute_progress_bars=True) as mpr:
        structure = mpr.get_structure_by_material_id(
            mpcode, conventional_unit_cell=is_conventional
        )

    return structure


def boxed_molecule(args):
    """Creates a boxed molecule from an input file"""
    from pymatgen.core import Molecule, Structure
    from pymatgen.io.vasp.inputs import Poscar

    # read the molecule
    box_molecule: Structure = Molecule.from_file(args.input).get_boxed_structure(
        a=args.vacuum, b=args.vacuum, c=args.vacuum, no_cross=args.no_cross
    )
    poscar = Poscar(box_molecule, sort_structure=True)

    if not args.output:
        print(poscar.get_str())
    else:
        poscar.write_file(f"{args.output}")


def convert_to_poscar(args):
    """Converts a file to a POSCAR file"""
    from pymatgen.io.vasp.inputs import Poscar

    structure = read(args.input)
    poscar = Poscar(structure, sort_structure=args.sort)

    if not args.output:
        print(poscar.get_str())
    else:
        poscar.write_file(f"{args.output}")

    return poscar

def list_selective_dynamics(args):
    """Lists atoms and their selective dynamics"""
    import pandas as pd
    from pymatgen.io.vasp.inputs import Poscar

    poscar = Poscar.from_file(args.input)

    atoms = read(args.input).get_chemical_symbols()
    dynamics = poscar.selective_dynamics
    heights = poscar.structure.cart_coords[:, 2]

    if dynamics is None:
        dynamics = [[True, True, True]] * len(atoms)

    atom_tuples = [
        (atom, index, sd, height)
        for atom, index, sd, height in zip(
            atoms, range(0, len(atoms)), dynamics, heights
        )
    ]
    # create a dataframe
    df = pd.DataFrame(
        atom_tuples, columns=["atom", "index", "selective_dynamics", "height"]
    )

    if not args.output:
        # print without index and print the whole dataframe
        with pd.option_context("display.max_rows", None, "display.max_columns", None):
            print(df.to_string(index=False))

    else:
        df.to_csv(args.output, index=False)



def calculate_rdf(coordinates, bins=1000, r_max=None):
    """
    Calculate the radial distribution function (RDF) from a set of XYZ coordinates.

    Parameters:
    - coordinates: numpy array of shape (n_particles, 3) representing XYZ coordinates.
    - bins: number of bins for histogram.
    - r_max: maximum distance to consider.

    Returns:
    - bin_centers: centers of the bins.
    - rdf: radial distribution function.
    """
    
    if r_max is None:
        r_max = np.max(np.linalg.norm(coordinates - coordinates[0], axis=1))

    # Calculate pairwise distances
    distances = np.sqrt(
        np.sum((coordinates[:, np.newaxis] - coordinates) ** 2, axis=-1)
    )

    # Create histogram
    hist, bin_edges = np.histogram(distances, bins=bins, range=(0, r_max))

    # Calculate bin centers and RDF
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    rdf = hist / (4 * np.pi * bin_centers**2 * (r_max / bins))

    # Remove the first bin because it is always 0
    bin_centers = bin_centers[1:]
    rdf = rdf[1:]

    # Normalize RDF
    rdf /= np.sum(rdf)

    return bin_centers, rdf


def plot_radial_distribution_function(args):
    """Plots the radial distribution function of a structure"""
    import plotly.graph_objects as go
    from pymatgen.core import Molecule

    coords = Molecule.from_file(args.input).cart_coords
    bin_centers, rdf = calculate_rdf(coords)

    # Plot RDF
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=bin_centers, y=rdf, mode="lines"))
    fig.update_layout(
        title="Radial distribution function",
        xaxis_title="Distance (Å)",
        yaxis_title="RDF",
        template="plotly_white",
    )

    if not args.output:
        fig.show()

    else:
        fig.write_image(args.output)


def compare_rdf(args):
    """Compares the radial distribution function of two structures"""
    import plotly.graph_objects as go
    from pymatgen.core import Molecule

    structures = [Molecule.from_file(file) for file in args.input]

    fig = go.Figure()
    for structure in structures:
        coords = structure.cart_coords
        bin_centers, rdf = calculate_rdf(coords)

        # Plot RDF
        fig.add_trace(
            go.Scatter(x=bin_centers, y=rdf, mode="lines", name=structure.formula)
        )

    fig.update_layout(
        title="Radial distribution function",
        xaxis_title="Distance (Å)",
        yaxis_title="RDF",
        template="plotly_white",
    )

    if not args.output:
        fig.show()

    else:
        fig.write_image(args.output)

def calculate_bond_angles(coords):
    """
    Calculate the bond angles from a set of XYZ coordinates.

    Parameters:
    - coords: numpy array of shape (n_particles, 3) representing XYZ coordinates.

    Returns:
    - bond_angles: bond angles.
    """
    from scipy.spatial import Delaunay

    # Calculate the Delaunay triangulation
    tri = Delaunay(coords)

    # Get the indices of the vertices of the tetrahedrons
    tetra = tri.simplices

    # Calculate the bond vectors
    bond_vectors = []
    for tet in tetra:
        for i in range(4):
            for j in range(i + 1, 4):
                bond_vectors.append(coords[tet[j]] - coords[tet[i]])

    bond_vectors = np.array(bond_vectors)

    # Calculate the bond angles
    bond_angles = []
    for i in range(0, len(bond_vectors), 3):
        a, b, c = bond_vectors[i : i + 3]
        bond_angles.append(
            np.degrees(
                np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
            )
        )
        bond_angles.append(
            np.degrees(
                np.arccos(np.dot(b, c) / (np.linalg.norm(b) * np.linalg.norm(c)))
            )
        )
        bond_angles.append(
            np.degrees(
                np.arccos(np.dot(c, a) / (np.linalg.norm(c) * np.linalg.norm(a)))
            )
        )

    return np.array(bond_angles)

def angular_distribution_function(bond_angles, bins=1000, r_max=2):
    """Calculate the angular distribution function (ADF) from a set of XYZ coordinates."""

    # Create histogram
    hist, bin_edges = np.histogram(bond_angles, bins=bins, range=(0, r_max))
    
    # Calculate bin centers and ADF
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    adf = hist / (4 * np.pi * bin_centers**2 * (r_max / bins))
    

    # Normalize ADF
    adf /= np.sum(adf)
    
    return bin_centers, adf



def compare_adf(args):
    """Compares the angular distribution function of two structures"""
    import plotly.graph_objects as go
    from pymatgen.core import Molecule

    structures = [Molecule.from_file(file) for file in args.input]

    fig = go.Figure()
    for structure in structures:
        coords = structure.cart_coords
        bond_angles = calculate_bond_angles(coords)
        bin_centers, adf = angular_distribution_function(bond_angles)

        # Plot ADF
        fig.add_trace(
            go.Scatter(x=bin_centers, y=adf, mode="lines", name=structure.formula)
        )

    fig.update_layout(
        title="Angular distribution function",
        xaxis_title="Angle (°)",
        yaxis_title="ADF",
        template="plotly_white",
    )

    if not args.output:
        fig.show()

    else:
        fig.write_image(args.output)
        
def structure_from_pubchem(cid: int, output: str = None):
    """Gets a structure from PubChem"""
    from pymatgen.core import Molecule
    import pubchempy as pcp
    compound = pcp.Compound.from_cid(cid, record_type='3d')
    
    #convert json data to xyz
    compound_dict = compound.to_dict(properties=['atoms'])
    symbols_and_coordinates = [(atom['element'], atom['x'], atom['y'], atom['z']) for atom in compound_dict['atoms']]

    #create a molecule
    species = [symbol for symbol, *_ in symbols_and_coordinates]
    coords = np.array([coord for _, *coord in symbols_and_coordinates])
    molecule = Molecule(species, coords)
    
    if not output:
        print(molecule.to(fmt="xyz"))
    else:
        molecule.to(fmt="xyz", filename=output)

def query_structure(args):
    """Queries XYZ file from PubChem"""
    #output format is {cid}_{output}
    if args.output:
        output_basenames = [f"{cid}_{args.output}" for cid in args.input]
        for cid, output in zip(args.input, output_basenames):
            structure_from_pubchem(cid, output)
    else:
        for cid in args.input:
            structure_from_pubchem(cid)
            


def run(args):
    functions = {
        "sort": sort_poscar,
        "rdf": compare_rdf,
        "adf": compare_adf,
        "query": query_structure,       
    }

    for arg, func in functions.items():
        if getattr(args, arg):
            func(args)
