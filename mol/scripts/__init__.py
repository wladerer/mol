import sys


def optimize(subparsers):
    subp_optimize = subparsers.add_parser('optimize', help="Optimize molecular geometry using UHF, RHF, or DFT")
    io_group = subp_optimize.add_argument_group("Input/Output")
    generic_group = subp_optimize.add_argument_group("Generic Information")
    method_group = subp_optimize.add_argument_group("Method-Specific Options")
    convergence_group = subp_optimize.add_argument_group("Convergence Options")

    # Input/Output arguments
    io_group.add_argument("input", default=sys.stdin, help="Input geometry")
    io_group.add_argument("-o", "--output", type=str, help="output filename")
    io_group.add_argument("--stdout", type=str, help="standard output stream")
    io_group.add_argument("--logfile", type=str, help="log filename")

    # Generic Information arguments
    generic_group.add_argument("-c", "--charge", type=int, help="charge")
    generic_group.add_argument("-s", "--spin", type=int, help="spin")
    generic_group.add_argument("-m", "--multiplicity", type=int, help="multiplicity")

    # Method-Specific Options arguments
    method_group.add_argument("-t", "--method", type=str, help="method", choices=["rhf", "uhf", "rohf", "ghf", "rks", "uks"])
    method_group.add_argument("-x", "--xc", type=str, help="functional")
    method_group.add_argument("-b", "--basis", type=str, help="basis")

    # Convergence Options arguments
    convergence_group.add_argument("-E", "--maxscf", type=int, help="maximum number of SCF steps", default=100)
    convergence_group.add_argument("-n", "--maxgeom", type=int, help="maximum number of ionic optimization steps", default=100)
    convergence_group.add_argument("-e", "--total-energy", type=float, help="convergence energy (eV)")
    convergence_group.add_argument("-g", "--grms", type=float, help="convergence grms (eV/Angstrom)")
    convergence_group.add_argument("-G", "--gmax", type=float, help="convergence gmax (eV/Angstrom)")
    convergence_group.add_argument("-d", "--drms", type=float, help="convergence drms (Angstrom)")
    convergence_group.add_argument("-D", "--dmax", type=float, help="convergence dmax (Angstrom)")


def structure(subparsers):
    subp_structure = subparsers.add_parser("structure", help="Generate and update structure files")

    subp_structure.add_argument(
        "input", default=sys.stdin, help="Input geometry or PubChem CID", nargs="+")
    subp_structure.add_argument(
        "--rdf",
        action="store_true",
        help="Plot the radial distribution function of a structure or multiple structures",
    )
    subp_structure.add_argument("--adf", action="store_true", help="Plot the angular distribution function of a structure or multiple structures")
    subp_structure.add_argument("-o", "--output", type=str, help="Output file name")
    subp_structure.add_argument("--sort", action="store_true", help="Sort atoms")
    subp_structure.add_argument(
        "-l", "--list", action="store_true", help="List the positions of the atoms"
    )
    subp_structure.add_argument("-q", "--query", action="store_true", help="Query a structure from PubChem")

    
def uvis(subparsers):
    
    subp_uvis = subparsers.add_parser("uvis", help="Plot UV-Vis spectra from TDDFT calculations")
    subp_uvis.add_argument("-o", "--output", type=str, help="Output file name")
    subp_uvis.add_argument("input", default=sys.stdin, help="Input geometry", nargs="+")
    subp_uvis.add_argument("-x", "--xc", type=str, help="Functional")
    subp_uvis.add_argument("-s", "--states", type=int, help="Number of states", default=15)
    subp_uvis.add_argument("-w", "--spectral-width", type=float, help="Spectral width", default=5)
    subp_uvis.add_argument("-k", "--key", type=str, help="Color key")
    
    #argument group for functional options
    functional_group = subp_uvis.add_argument_group("Functional Options")
    functional_group.add_argument("--plot", action="store_true", help="Plot the UV-Vis spectra")
    functional_group.add_argument("--save", action="store_true", help="Save the UV-Vis spectra to a csv file")
    functional_group.add_argument("--read", action="store_true", help="Plot the UV-Vis spectra from a csv file")
    

def ir(subparsers):
    
    subp_ir = subparsers.add_parser("ir", help="Plot IR spectra from KSDFT calculations")
    subp_ir.add_argument("-o", "--output", type=str, help="Output file name")
    subp_ir.add_argument("input", default=sys.stdin, help="Input geometry", nargs="+")
    subp_ir.add_argument("-k", "--key", type=str, help="Color key")
    subp_ir.add_argument("-S", "--scale", type=float, help="Scale factor", default=1.0)
    subp_ir.add_argument("-w", "--fwhw", type=float, help="Full width at half width", default=100.0)
    
    #argument group for functional options
    functional_group = subp_ir.add_argument_group("Functional Options")
    functional_group.add_argument("--plot", action="store_true", help="Plot the IR spectra")
    functional_group.add_argument("--save", action="store_true", help="Save the IR spectra to a csv file")
    functional_group.add_argument("--read", action="store_true", help="Plot the IR spectra from a csv file")
    
    # Generic Information arguments
    generic_group = subp_ir.add_argument_group("Generic Information")
    generic_group.add_argument("-c", "--charge", type=int, help="charge")
    generic_group.add_argument("-s", "--spin", type=int, help="spin")
    generic_group.add_argument("-m", "--multiplicity", type=int, help="multiplicity")
    
    # Method-Specific Options arguments
    method_group = subp_ir.add_argument_group("Method-Specific Options")
    method_group.add_argument("-t", "--method", type=str, help="method", choices=["rks", "uks"])
    method_group.add_argument("-x", "--xc", type=str, help="functional")
    method_group.add_argument("-b", "--basis", type=str, help="basis")

def setup(subparsers):
    scripts = [
        structure,
        optimize,
        uvis,
        ir
    ]
    for script in scripts:
        script(subparsers)
