from mcp.server.fastmcp import FastMCP
from mcp.types import TextContent, ImageContent, BlobResourceContents
import logging

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

# --- NEW TOOL FOR PLOTTING ---
@mcp.tool(
    name="plot_energy_scan_image",
    description="Generates a plot of energy versus bond length from provided data and returns it as a Base64 encoded PNG image.",
    parameters=[
        ToolParameter(name="bond_lengths", type="array", items={"type": "number"}, description="List of bond lengths (x-axis data).", required=True),
        ToolParameter(name="energies", type="array", items={"type": "number"}, description="List of corresponding energies (y-axis data).", required=True),
        ToolParameter(name="title", type="string", description="Title of the plot.", required=False, default="Energy vs. Bond Length Scan"),
        ToolParameter(name="xlabel", type="string", description="Label for the x-axis.", required=False, default="Bond Length (Angstroms)"),
        ToolParameter(name="ylabel", type="string", description="Label for the y-axis.", required=False, default="Energy (Hartree)")
    ],
    returns=ToolParameter(
        type="object",
        properties={
            "image_base64": ToolParameter(type="string", description="Base64 encoded PNG image data."),
            "format": ToolParameter(type="string", description="Format of the image (e.g., 'png').")
        },
        description="A dictionary containing the Base64 encoded image and its format."
    )
)
def plot_energy_scan_image_mcp(bond_lengths: List[float], energies: List[float],
                                title: str = "Energy vs. Bond Length Scan",
                                xlabel: str = "Bond Length (Angstroms)",
                                ylabel: str = "Energy (Hartree)") -> Dict[str, str]:
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
