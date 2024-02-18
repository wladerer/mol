from pyscf import gto
from pyscf.prop import infrared
from pyscf.hessian import thermo

def plot_ir_from_mf(mf, mol: gto.Mole, temp: float = 298.15, pressure: float = 101325):
    """
    Plots the infrared spectrum from a mean-field calculation.

    Args:
        mf (object): The mean-field calculation object.
        mol (gto.Mole): The molecular system.
        temp (float, optional): The temperature in Kelvin. Defaults to 298.15.
        pressure (float, optional): The pressure in Pascal. Defaults to 101325.
    """
    mf_ir = infrared.rks.Infrared(mf).run()
    mf_ir.summary()
    
    freq_au = mf_ir.vib_dict["freq_au"]
    thermo_results = thermo.thermo(mf, freq_au, temperature=temp, pressure=pressure)
    thermo.dump_thermo(mol, thermo_results)
    
    fig, ax, ax2 = mf_ir.plot_ir(w=100, scale=0.956)
    ax.set_title(r"Infrared Spectrum")
    fig.tight_layout()
    fig.show()

