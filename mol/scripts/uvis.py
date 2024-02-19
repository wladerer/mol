from mol.utils.spectra import plot_uv_vis_from_df
from mol.utils.io import geometry_from_file
from pyscf import gto, dft, tddft
import pandas as pd
import numpy as np
from scipy.constants import physical_constants

ha_2_ev = 1/physical_constants["electron volt-hartree relationship"][0]

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def run_spectral_analysis(mol: gto.Mole, xc: str, states: int = 15, spectral_width: float = 0.1) -> pd.DataFrame:

    # Ground State DFT
    mf = dft.RKS(mol, xc=xc).run()

    # Excited State DFT
    mytd = tddft.TDDFT(mf)
    mytd.nstates = states
    mytd.max_space = 100
    mytd.max_cycle = 200
    mytd.kernel()
    mytd.analyze()
    osc_strengths = mytd.oscillator_strength()[:states-5]

    # Convolve lineshapes to make spectra
    energies_ev = mytd.e[:states-5]*ha_2_ev
    x_range = np.linspace(energies_ev.min()*0.9, energies_ev.max()*1.1, num=1000)
    intensity = np.zeros(x_range.size)

    for e, f in zip(energies_ev, osc_strengths):
        intensity += gaussian(x_range, e, spectral_width) * f

    # Rough Normalization
    dx = (x_range[-1] - x_range[0])/x_range.size
    area = (intensity*dx).sum()
    intensity /= area

    # x="Excitation Energy (eV)", y="Intensity",
    df = pd.DataFrame({"Excitation Energy (eV)": x_range, "Intensity": intensity})

    return df

def run_and_plot(args):
    
    mols = [ geometry_from_file(file) for file in args.input ]
    df = pd.concat([ run_spectral_analysis(mol, args.xc, args.states, args.spectral_width) for mol in mols ])
    # formulae = [] for later, if we want to add the formulae to the plot as labels
    plot_uv_vis_from_df(df, color_key=args.key)
    
    
def run_and_save(args):
    
    mols = [ geometry_from_file(file) for file in args.input ]
    df = pd.concat([ run_spectral_analysis(mol, args.xc, args.states, args.spectral_width) for mol in mols ])
    df.to_csv(args.output, index=False)
    
def plot_from_csv(args):
    
    df = pd.read_csv(args.input)
    plot_uv_vis_from_df(df, color_key=args.key)
    

def run(args):
    
    functions = {
        "plot": run_and_plot,
        "save": run_and_save,
        "read": plot_from_csv
    }
    
    for key, value in functions.items():
        if getattr(args, key):
            value(args)
            break
    