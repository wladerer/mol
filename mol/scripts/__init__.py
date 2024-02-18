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
    convergence_group.add_argument("-n", "--maxsteps", type=int, help="maximum number of optimization steps")
    convergence_group.add_argument("-e", "--total-energy", type=float, help="convergence energy (eV)")
    convergence_group.add_argument("-g", "--grms", type=float, help="convergence grms (eV/Angstrom)")
    convergence_group.add_argument("-G", "--gmax", type=float, help="convergence gmax (eV/Angstrom)")
    convergence_group.add_argument("-d", "--drms", type=float, help="convergence drms (Angstrom)")
    convergence_group.add_argument("-D", "--dmax", type=float, help="convergence dmax (Angstrom)")


def structure(subparsers):
    subp_structure = subparsers.add_parser("structure", help="Generate and update structure files")

    subp_structure.add_argument(
        "input", default=sys.stdin, help="Input geometry")
    subp_structure.add_argument(
        "--rdf",
        action="store_true",
        help="Plot the radial distribution function of a structure",
    )
    subp_structure.add_argument("-o", "--output", type=str, help="Output file name")
    subp_structure.add_argument("--sort", action="store_true", help="Sort atoms")
    subp_structure.add_argument(
        "-l", "--list", action="store_true", help="List the positions of the atoms"
    )
    subp_structure.add_argument(
        "--compare",
        nargs="+",
        type=str,
        help="Compare the radial distribution of multiple structures",
    )


def setup(subparsers):
    scripts = [
        structure,
        optimize,
    ]
    for script in scripts:
        script(subparsers)
