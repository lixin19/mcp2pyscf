from mcp.server.fastmcp import FastMCP
from mcp.types import TextContent, ImageContent, BlobResourceContents
import logging
from typing import Literal, List, Dict, Union # Import new types
from pyscf_tools import smiles_to_xyz, run_bond_stretch_scan
from plot import plot_energy_scan # Import new functions


# Set up logging (this just prints messages to your terminal for debugging)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create the MCP server object
mcp = FastMCP()

@mcp.tool()
def pyscf_rhf_energy(atom, basis):
    """
    Run a restricted Hartree-Fock calculation (RHF) with PySCF
    
    input:
        atom: the geometry of the molecule in Angstrom
    """
    from pyscf import gto, scf

    mol = gto.M(
        #atom = 'H 0 0 0; H 0 0 1.1',  # in Angstrom
        atom=atom,
        basis = basis,
        symmetry = True,
    )
    mf = scf.HF(mol)
    mf.kernel()
    e_tot = mf.e_tot
    #return e_tot
    return TextContent(type="text", text=str(e_tot))

@mcp.tool()
def scan_pes_rhf():
    """
    Scan a potential energy curve with RHF
    """
    import numpy as np
    from pyscf import scf
    from pyscf import gto

    ehf = []

    # hard code basis for now
    basis = "sto-3g"
        
    for b in np.arange(0.7, 1.51, 0.1):
        # hard code molecule for now
        atom = f"H 0 0 0; H 0 0 {b}"  # in Angstrom

        result = pyscf_rhf_energy(atom,basis)
        ehf.append(float(result.text))  # convert string back to float

    print('R     E(HF)')
    for i, b in enumerate(np.arange(0.7, 1.51, 0.1)):
        print('%.2f  %.8f' % (b, ehf[i]))

    return TextContent(type="text", text=str(ehf))

@mcp.tool()
def generate_pyscf_geom_input(smiles_string: str) -> str:
    """
    Generates a geometry string in PySCF format for a PySCF calculation from a SMILES string.

    Args:
        smiles_string: The SMILES string of the molecule (e.g., "CCO" for ethanol).

    Returns:
        A string containing the Python script for PySCF.
    """
    try:
        return smiles_to_xyz(smiles_string)
    except ValueError as e:
        return f"Error: {e}. Please provide a valid SMILES string."

# --- NEW TOOL FOR BOND STRETCH SCAN ---
@mcp.tool()
def run_bond_stretch_calculation_mcp(smiles_string: str, atom1_idx: int, atom2_idx: int,
                                    start_dist: float, end_dist: float, num_points: int, basis: str = "sto-3g") -> Dict[str, List[float]]:
    """
    Performs a series of PySCF Hartree-Fock/STO-3G energy calculations
    for varying a specified bond length in a molecule and returns bond lengths and energies.

    Args:
        smiles_string: SMILES string of the molecule.
        atom1_idx: 0-indexed ID of the first atom in the bond.
        atom2_idx: 0-indexed ID of the second atom in the bond.
        start_dist: Starting bond distance in Angstroms.
        end_dist: Ending bond distance in Angstroms.
        num_points: Number of points to calculate in the scan.

    Returns:
        A dictionary containing 'bond_lengths' (list of floats) and 'energies' (list of floats).
    """
    try:
        # Call your actual PySCF function
        results = run_bond_stretch_scan(smiles_string, atom1_idx, atom2_idx, start_dist, end_dist, num_points, basis)
        return results
    except Exception as e:
        return {"error": str(e), "bond_lengths": [], "energies": []}


# --- NEW TOOL FOR PLOTTING ---
@mcp.tool()
def plot_energy_scan_image_mcp(bond_lengths: List[float], energies: List[float],
                                title: str = "Energy vs. Bond Length Scan",
                                xlabel: str = "Bond Length (Angstroms)",
                                ylabel: str = "Energy (Hartree)") -> Dict[str, str]:
    """
    Generates a plot of energy versus bond length from provided data and returns it as a Base64 encoded PNG image.

    Args:
        bond_lengths: List of bond lengths (x-axis data).
        energies: List of energies (y-axis data).
        title: Title of the plot.
        xlabel: Label for the x-axis.
        ylabel: Label for the y-axis.

    Returns:
        A dictionary containing 'image_base64' (Base64 encoded PNG image data) and 'format' (e.g., 'png').
    """
    try:
        # Call your actual plotting function
        img_b64 = plot_energy_scan(bond_lengths, energies, title, xlabel, ylabel)
        return {"image_base64": img_b64, "format": "png"}
    except Exception as e:
        return {"error": str(e), "image_base64": "", "format": ""}



# This is the main entry point for your server
def main():
    logger.info('Starting your-new-server')
    mcp.run('stdio')

if __name__ == "__main__":
    main()
