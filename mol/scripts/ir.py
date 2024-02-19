from pyscf.prop import infrared
from pyscf import dft
from pyscf.hessian import thermo

from mol.utils.io import geometry_from_file

import matplotlib.pyplot as plt


def plot_ir_from_mol(mol, mf: dft.RKS) -> plt.Figure:
    mf_ir = infrared.rks.Infrared(mf).run()
    mf_ir.summary()

    thermo.dump_thermo(mol, thermo.thermo(mf, mf_ir.vib_dict["freq_au"], 298.15, 101325))

    fig, ax, ax2 = mf_ir.plot_ir(w=100, scale=0.956)
    ax.set_title(r"Infrared Spectrum")
    fig.tight_layout()
    
    return fig

def run_and_plot(args):
    
    mols = [ geometry_from_file(file) for file in args.input ]
    mfs = [ dft.RKS(mol).run() for mol in mols]
    figs = [ plot_ir_from_mol(mol, mf) for mol, mf in zip(mols, mfs) ]
    
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
    }
    
    for key, value in functions.items():
        if getattr(args, key):
            value(args)
            break
