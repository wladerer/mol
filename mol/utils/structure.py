import numpy as np
import pandas as pd
from ase.io import read
from scipy.spatial import distance_matrix


def xyz_to_dataframe(file: str):
    """Reads an xyz file and returns a pandas dataframe."""
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        lines = [line.split() for line in lines]
        lines = [line for line in lines if len(line) == 4]
        df = pd.DataFrame(lines, columns=["atom", "x", "y", "z"])
        df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(float)

        return df


def sort_df_by_height(df: pd.DataFrame):
    """Sorts a dataframe by the z coordinate."""
    df = df.sort_values(by=["z"], ascending=False)
    df = df.reset_index(drop=True)

    return df



def df_to_distance_matrix(df: pd.DataFrame, n: int = 12):
    """Reads a structure file and returns a distance matrix"""
    data = sort_df_by_height(df)
    dist = distance_matrix(data.iloc[:n, 1:], data.iloc[:n, 1:])
    dist = pd.DataFrame(dist)
    atom_labels = data.iloc[:n, 0].values
    dist.index = atom_labels
    dist.columns = atom_labels

    return dist


def calculate_rdf(coordinates, bins=1000, r_max=15):
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


def plot_radial_distribution_function(
    input_files: list[str], labels: list[str], output_file: str = None
):
    """Plots the radial distribution function of a structure"""
    import plotly.graph_objects as go
    from pymatgen.core import Molecule

    fig = go.Figure()

    for input_file, label in zip(input_files, labels):
        coords = Molecule.from_file(input_file).coords
        bin_centers, rdf = calculate_rdf(coords)

        # Add trace for each input file
        fig.add_trace(go.Scatter(x=bin_centers, y=rdf, mode="lines", name=label))

    fig.update_layout(
        title="Radial distribution function",
        xaxis_title="Distance (Å)",
        yaxis_title="RDF",
        template="plotly_white",
    )

    if not output_file:
        fig.show()
    else:
        fig.write_image(output_file)


def angular_distribution_function(coordinate, bins=1000, r_max=180):
    """Calculate the angular distribution function (ADF) from a set of XYZ coordinates."""
    # Calculate the bond angles between all atoms
    bond_angles = []
    for i in range(len(coordinate)):
        for j in range(i + 1, len(coordinate)):
            for k in range(j + 1, len(coordinate)):
                v1 = coordinate[i] - coordinate[j]
                v2 = coordinate[k] - coordinate[j]
                angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
                bond_angles.append(np.degrees(angle))
                
    # Create histogram
    hist, bin_edges = np.histogram(bond_angles, bins=bins, range=(0, r_max))
    
    # Calculate bin centers and ADF
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    adf = hist / (4 * np.pi * bin_centers**2 * (r_max / bins))
    
    # Remove the first bin because it is always 0
    bin_centers = bin_centers[1:]
    adf = adf[1:]
    
    # Normalize ADF
    adf /= np.sum(adf)
    
    return bin_centers, adf

def plot_angular_distribution_function(input_files: list[str], labels: list[str], output_file: str = None):
    """Plots the angular distribution function of a list of structures"""
    
    import plotly.graph_objects as go
    from pymatgen.core import Molecule
    
    fig = go.Figure()
    
    for input_file, label in zip(input_files, labels):
        coords = Molecule.from_file(input_file).coords
        bin_centers, adf = angular_distribution_function(coords)
        
        # Add trace for each input file
        fig.add_trace(go.Scatter(x=bin_centers, y=adf, mode="lines", name=label))
        
    fig.update_layout(
        title="Angular distribution function",
        xaxis_title="Angle (°)",
        yaxis_title="ADF",
        template="plotly_white"
    )
    
    if not output_file:
        fig.show()
    else:
        fig.write_image(output_file)


### symmetry tool section


def get_symmetry_operations(file: str):
    """Returns the symmetry operations of a POSCAR or CONTCAR file."""
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointgroupAnalyzer

    structure = Molecule.from_file(file)
    sga = PointgroupAnalyzer(structure)
    symmetry_operations = sga.get_symmetry_operations()

    return symmetry_operations


def plot_atomic_drift(initial_file: str, final_file: str, output_file: str = None):
    """Plots the atomic drift between two structures"""
    import sisl
    import sisl.viz

    intial_atoms = read(initial_file)
    final_atoms = read(final_file)
    drift = final_atoms.get_positions() - intial_atoms.get_positions()
    geom = sisl.get_sile(initial_file).read_geometry()
    plot = geom.plot()

    plot.update_inputs(
        arrows={"data": drift, "name": "Drift", "color": "orange", "width": 2},
        axes="xyz",
    )

    plot.show()
