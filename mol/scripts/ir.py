from pyscf.prop import infrared
from pyscf.prop.infrared.rhf import ir_point
from pyscf import dft, gto
from pyscf.hessian import thermo

from mol.utils.io import geometry_from_file

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def run_ir(mol: gto.Mole, mf: dft.UKS) -> infrared.uks.Infrared:
    mf_ir = infrared.uks.Infrared(mf).run()
    mf_ir.summary()

    thermo.dump_thermo(mol, thermo.thermo(mf, mf_ir.vib_dict["freq_au"], 298.15, 101325))
    
    return mf_ir



def get_plot_from_df(df: pd.DataFrame, fwhw: float = 100, scale: float = 1.0) -> plt.Figure:
    '''Creates a plot from pandas dataframe'''

    freq = df["freq_wavenumber"]
    ir_inten = df["ir_inten"]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.grid()
    
    x = np.linspace(0, np.max([4000., np.max(freq) + 5*fwhw]), 4000)
    ax.plot(x, [ir_point(xi, fwhw, freq, ir_inten) for xi in x], label="Molar Absorption Coefficient")
    ax.set_ylabel("Molar Absorption Coefficient (L mol$^{-1}$ cm$^{-1}$)")
    ax.set_xlabel("Vibration Wavenumber (cm$^{-1}$)")
    ax.legend(loc="upper left")

    ax2 = ax.twinx()
    for i in range(ir_inten.size):
        if i == 0:
            ax2.plot([freq[i], freq[i]], [0, ir_inten[i]], c="C2", linewidth=1, label="IR Intensity")
        else:
            ax2.plot([freq[i], freq[i]], [0, ir_inten[i]], c="C2", linewidth=1)
    ax2.set_ylabel("IR Intensity (km mol$^{-1}$)")
    ax2.legend(loc="upper right")
    fig.tight_layout()
    
    return fig, ax, ax2



def plot_ir(mf_ir: infrared.uks.Infrared, fwhw: float = 100, scale: float = 1.0) -> plt.Figure:

    fig, ax, ax2 = mf_ir.plot_ir(w=fwhw, scale=scale)
    ax.set_title(r"Infrared Spectrum")
    fig.tight_layout()
    
    return fig

def ir_to_df(mf_ir: infrared.uks.Infrared, scale: float, series: int) -> pd.DataFrame:
    
    freq = mf_ir.vib_dict["freq_wavenumber"] * scale
    freq = np.real(freq)
    
    ir_inten = mf_ir.ir_inten.copy()
    ir_inten[np.abs(np.imag(freq)) > 1e-10] = 0
    
    return pd.DataFrame({"freq_wavenumber": freq, "ir_inten": ir_inten, "series": series})

  
    
def run_and_save(args):
        
        mols = [ geometry_from_file(file) for file in args.input ]
        mfs = [ dft.UKS(mol).run() for mol in mols]
        mf_irs = [ run_ir(mol, mf) for mol, mf in zip(mols, mfs) ]
        dfs = [ ir_to_df(mf_ir, args.scale, i) for i, mf_ir in enumerate(mf_irs) ]
        
        
        if not args.output:
            for df in dfs:
                print(df)
        
        else:
            for i, df in enumerate(dfs):
                if len(dfs) > 1:
                    df.to_csv(f"{i+1}_{args.output}")
                else:
                    dfs[0].to_csv(args.output)


def read_and_plot(args):
    """Reads a cvs file and plots the IR spectra"""
    dfs = [ pd.read_csv(file) for file in args.input ]
    figs = [ get_plot_from_df(df, args.fwhw, args.scale)[0] for df in dfs ]
    
    if not args.output:
        for fig in figs:
            fig.show()
    else:
        for i, fig in enumerate(figs):
            if len(figs) > 1:
                for i, fig in enumerate(figs):
                    fig.savefig(f"{i+1}_{args.output}")
            else:
                figs[0].savefig(args.output)


def run_and_plot(args):
    
    mols = [ geometry_from_file(file) for file in args.input ]
    mfs = [ dft.UKS(mol).run() for mol in mols]
    mf_irs = [ run_ir(mol, mf) for mol, mf in zip(mols, mfs) ]
    figs = [ plot_ir(mf_ir, args.fwhw, args.scale) for mf_ir in mf_irs ]
    
    if not args.output:
        for fig in figs:
            fig.show()
    else:
        for i, fig in enumerate(figs):
            if len(figs) > 1:
                for i, fig in enumerate(figs):
                    fig.savefig(f"{i+1}_{args.output}")
            else:
                figs[0].savefig(args.output)

def run(args):
    
    functions = {
        "plot": run_and_plot,
        "save": run_and_save,
        "read": read_and_plot
    }
    
    for key, value in functions.items():
        if getattr(args, key):
            value(args)
            break
