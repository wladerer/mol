from pyscf.prop import infrared
from pyscf import dft, gto
from pyscf.hessian import thermo

from mol.utils.io import geometry_from_file

import matplotlib.pyplot as plt
import pandas as pd

def run_ir(mol: gto.Mole, mf: dft.UKS) -> infrared.uks.Infrared:
    mf_ir = infrared.uks.Infrared(mf).run()
    mf_ir.summary()

    thermo.dump_thermo(mol, thermo.thermo(mf, mf_ir.vib_dict["freq_au"], 298.15, 101325))
    
    return mf_ir

def plot_ir(mf_ir: infrared.uks.Infrared) -> plt.Figure:

    fig, ax, ax2 = mf_ir.plot_ir(w=100, scale=0.956)
    ax.set_title(r"Infrared Spectrum")
    fig.tight_layout()
    
    return fig

def save_ir(mf_ir: infrared.uks.Infrared, output: str):
    vib_dict = mf_ir.vib_dict
    vib_df = pd.DataFrame(vib_dict)
    vib_df.to_csv(output, index=False)
    
def run_and_save(args):
        
        mols = [ geometry_from_file(file) for file in args.input ]
        mfs = [ dft.UKS(mol).run() for mol in mols]
        mf_irs = [ run_ir(mol, mf) for mol, mf in zip(mols, mfs) ]
        [ save_ir(mf_ir, f"{i+1}_{args.output}") for i, mf_ir in enumerate(mf_irs) ]


def run_and_plot(args):
    
    mols = [ geometry_from_file(file) for file in args.input ]
    mfs = [ dft.UKS(mol).run() for mol in mols]
    mf_irs = [ run_ir(mol, mf) for mol, mf in zip(mols, mfs) ]
    figs = [ plot_ir(mf_ir) for mf_ir in mf_irs ]
    
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
        "save": run_and_save
    }
    
    for key, value in functions.items():
        if getattr(args, key):
            value(args)
            break
