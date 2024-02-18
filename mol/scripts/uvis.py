from mol.utils.irplot import plot_uv_vis_from_df
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


mol_2Li = geometry_from_file('/home/wladerer/research/F4TCNQ/cyano/opt/opt_2Li_cyanoF4TCNQ.xyz')
mol_2Na = geometry_from_file('/home/wladerer/research/F4TCNQ/cyano/opt/opt_2Na_cyanoF4TCNQ.xyz')
mol_2Rb = geometry_from_file('/home/wladerer/research/F4TCNQ/cyano/opt/opt_2Rb_cyanoF4TCNQ.xyz')

import concurrent.futures

def run_spectral_analysis_parallel(mol, xc, states, spectral_width):
    return run_spectral_analysis(mol, xc, states, spectral_width)

with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = []
    futures.append(executor.submit(run_spectral_analysis_parallel, mol_2Li, 'b3lyp', 12, 0.1))
    futures.append(executor.submit(run_spectral_analysis_parallel, mol_2Na, 'b3lyp', 12, 0.1))
    futures.append(executor.submit(run_spectral_analysis_parallel, mol_2Rb, 'b3lyp', 12, 0.1))

    df_2Li = futures[0].result()
    df_2Na = futures[1].result()
    df_2Rb = futures[2].result()

#combine all the dataframes
df_2Li['System'] = '2Li'
df_2Na['System'] = '2Na'
df_2Rb['System'] = '2Rb'

#add functionals to the dataframes

df_2Li['Functional'] = 'B3LYP'
df_2Na['Functional'] = 'B3LYP'
df_2Rb['Functional'] = 'B3LYP'

df = pd.concat([df_2Li, df_2Na, df_2Rb])
#save the dataframe to csv
df.to_csv('spectra.csv', index=False)

plot_uv_vis_from_df(df, color_key='System')

# run 2Ag

mol_2Ag = geometry_from_file('/home/wladerer/research/F4TCNQ/cyano/opt/opt_2Ag_cyanoF4TCNQ.xyz')

df_2Ag = run_spectral_analysis(mol_2Ag, 'b3lyp', 12, 0.1)
df_2Ag['System'] = '2Ag'
df_2Ag['Functional'] = 'B3LYP'
df_2Ag.to_csv('2Ag_spectra.csv', index=False)