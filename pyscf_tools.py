### python 3
# To generate sample pyscf input
#   can read in sample input from /sample/ directory
#   then add in geometry information
###

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, scf
from pyscf.lib import parameters
from plot import plot_energy_scan
from pyscf.geomopt.geometric_solver import optimize
import webbrowser

def smiles_to_xyz(smiles_string: str) -> str:
    """
    Generates geometry information in PySCF atom string format for the geometry for PySCF energy calculation
    from a SMILES string. Includes 3D geometry generation using RDKit.

    Args:
        smiles_string: The SMILES string of the molecule (e.g., "CCO" for ethanol).

    Returns:
        A string containing the Python script for PySCF.
    """
    # 1. Generate 3D geometry from SMILES using RDKit
    mol_rdkit = Chem.MolFromSmiles(smiles_string)
    if mol_rdkit is None:
        raise ValueError(f"Invalid SMILES string: {smiles_string}")

    mol_rdkit = Chem.AddHs(mol_rdkit) # Add explicit hydrogens
    AllChem.EmbedMolecule(mol_rdkit, AllChem.ETKDG()) # Generate 3D conformer
    AllChem.UFFOptimizeMolecule(mol_rdkit) # Optimize geometry with UFF (optional, but good for stability)

    # Convert RDKit molecule to PySCF atom string format
    pyscf_atom_string = ""
    conf = mol_rdkit.GetConformer()
    for i, atom in enumerate(mol_rdkit.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        pyscf_atom_string += f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f};\n"
    pyscf_atom_string = pyscf_atom_string.strip() # Remove trailing semicolon and newline

    return pyscf_atom_string


# --- FUNCTION FOR BOND STRETCHING ---
def run_bond_stretch_scan(smiles_string: str, atom1_idx: int, atom2_idx: int,
                          start_dist: float, end_dist: float, num_points: int, basis: str = "sto-3g") -> dict:
    """
    Performs a bond stretching scan (energy vs. bond length) for a specified bond
    using PySCF with HF/STO-3G.

    Args:
        smiles_string: SMILES string of the molecule.
        atom1_idx: 0-indexed integer for the first atom in the bond.
        atom2_idx: 0-indexed integer for the second atom in the bond.
        start_dist: Starting bond distance in Angstroms.
        end_dist: Ending bond distance in Angstroms.
        num_points: Number of points in the scan.

    Returns:
        A dictionary containing 'bond_lengths' (list of floats) and 'energies' (list of floats).
    """
    mol_rdkit = Chem.MolFromSmiles(smiles_string)
    if mol_rdkit is None:
        raise ValueError(f"Invalid SMILES string: {smiles_string}")
    mol_rdkit = Chem.AddHs(mol_rdkit)
    AllChem.EmbedMolecule(mol_rdkit, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol_rdkit) # Optional: initial optimization

    bond_lengths = np.linspace(start_dist, end_dist, num_points).tolist()
    energies = []

    conf = mol_rdkit.GetConformer()

    for dist in bond_lengths:
        # Get current positions
        coords = conf.GetPositions()
        
        # Calculate current bond vector and its length
        vec = coords[atom2_idx] - coords[atom1_idx]
        current_dist = np.linalg.norm(vec)
        
        # Scale the vector to the desired new distance
        if current_dist == 0: # Avoid division by zero if atoms are coincident
             raise ValueError("Atoms are coincident, cannot perform bond stretch.")
        
        scale_factor = dist / current_dist
        
        # Apply scaling to move atom2_idx relative to atom1_idx
        # A more robust approach for geometry manipulation would involve internal coordinates
        # or a dedicated geometry optimization library like geomeTRIC.
        # For a hackathon, a simple scaling is often sufficient for a demo.
        new_vec = vec * scale_factor
        
        # Update atom2_idx position, keeping atom1_idx fixed for simplicity
        coords[atom2_idx] = coords[atom1_idx] + new_vec
        
        # Construct PySCF molecule object
        pyscf_atom_list = []
        for i, atom in enumerate(mol_rdkit.GetAtoms()):
            pyscf_atom_list.append([atom.GetSymbol(), coords[i][0], coords[i][1], coords[i][2]])

        mol = gto.M(atom=pyscf_atom_list, basis='sto-3g')

        # Run HF calculation
        mf = scf.RHF(mol).run(verbose=0) # verbose=0 to suppress extensive output

        energies.append(mf.e_tot)

    return {"bond_lengths": bond_lengths, "energies": energies}

# --- FUNCTION FOR GEOMETRY OPTIMIZATION ---
def optimize_molecule(smiles_string: str, maxsteps: int = 100) -> str:
    """
    Performs a geometry optimization for a molecule using PySCF with HF/STO-3G.

    Args:
        smiles_string: The SMILES string of the molecule.

    Returns:
        A dictionary containing:
        - 'optimized_energy': The final optimized energy in Hartree.
        - 'optimized_geometry_xyz': The optimized geometry as an XYZ string.
    """
    # 1. Generate initial 3D geometry from SMILES using RDKit
    mol_rdkit = Chem.MolFromSmiles(smiles_string)
    if mol_rdkit is None:
        raise ValueError(f"Invalid SMILES string: {smiles_string}")
    mol_rdkit = Chem.AddHs(mol_rdkit)
    AllChem.EmbedMolecule(mol_rdkit, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol_rdkit) # Initial classical optimization for a reasonable start

    # Convert RDKit molecule to PySCF atom list format
    # PySCF's optimize expects atom as a list of lists or an array
    initial_atom_list = []
    conf = mol_rdkit.GetConformer()
    for i, atom in enumerate(mol_rdkit.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        initial_atom_list.append([atom.GetSymbol(), pos.x, pos.y, pos.z])

    mol_pyscf = gto.M(
        atom = initial_atom_list,
        basis = 'sto-3g'
    )

    # 2. Perform geometry optimization
    # This creates a Gradients object (Gradients from an SCF object)
    mf = scf.RHF(mol_pyscf).run(verbose=0) # Run initial SCF calculation
    
    # Use PySCF's built-in optimizer (Newton-Raphson by default)
    # The .kernel() method of the Gradients object performs the optimization
    optimized_mol = optimize(mf, maxsteps) # Pass the mean-field object (mf) to geom.optimize

    #optimized_energy = optimized_mol.e_tot # The optimized molecule object has the final energy

    # Convert optimized PySCF molecule geometry to XYZ string
    optimized_geometry_xyz = f"{optimized_mol.natm}\n"
    optimized_geometry_xyz += "Optimized geometry (HF/STO-3G)\n"
    for i in range(optimized_mol.natm):
        sym = optimized_mol.atom_symbol(i)
        xyz_bohr = optimized_mol.atom_coord(i)
        xyz = xyz_bohr * parameters.BOHR # units.BOHR is 0.529177... Angstroms/Bohr
        optimized_geometry_xyz += f"{sym} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n"
    print(f"Optimized geometry for {smiles_string}:\n{optimized_geometry_xyz}")
    return optimized_geometry_xyz



def generate_3d_viewer_html(xyz_string: str, title: str = "Molecular Visualization") -> str:
    """
    Generates an HTML string with an embedded 3Dmol.js viewer for a given XYZ molecular string.

    Args:
        xyz_string: The molecular geometry in XYZ format.
        title: Title for the HTML page.

    Returns:
        A string containing the complete HTML for the 3D viewer.
    """
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; overflow: hidden; }}
        #viewer_container {{ width: 100vw; height: 100vh; position: relative; }}
    </style>
</head>
<body>
    <div id="viewer_container"></div>
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            let element = document.getElementById('viewer_container');
            let viewer = $3Dmol.createViewer( element, {{ backgroundColor: 'white' }} );

            const xyz_data = `{xyz_string}`; // Use backticks for multi-line string literal

            viewer.addModel( xyz_data, "xyz" );
            // --- MODIFIED STYLE FOR BONDS AND ATOMS ---
            viewer.setStyle({{ 
                'stick': {{ 'radius': 0.2, 'color': 'lightgray' }}, // Increased radius for thicker bonds, added color
                'sphere': {{ 'scale': 0.30, 'colorscheme': 'jmol' }} // Slightly larger atoms, using Jmol color scheme
            }});
            // You can also try 'line' for a simpler representation:
            // viewer.setStyle({{ 'line': {{ 'linewidth': 3 }} }});
            // Or 'bond' for explicit bond rendering (often used with explicit bond info, not just XYZ)
            // viewer.setStyle({{ 'bond': {{ 'color': 'black', 'radius': 0.2 }} }});
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
</body>
</html>
"""
    return html_content


# Example usage for the new functions (for local testing)
if __name__ == "__main__":
    # # Test bond stretching
    # try:
    #     print("\n--- Testing Bond Stretch Scan and Plotting ---")
    #     scan_results = run_bond_stretch_scan(
    #         smiles_string="O", # Water
    #         atom1_idx=0, atom2_idx=1, # O-H bond. RDKit usually has O=0, H1=1, H2=2
    #         start_dist=0.8, end_dist=1.5, num_points=10, basis='sto-3g'
    #     )
    #     print("\nBond Scan Results (first few):")
    #     print(f"  Lengths: {scan_results['bond_lengths'][:3]}...")
    #     print(f"  Energies: {scan_results['energies'][:3]}...")

    #     # Test plotting
    #     img_b64 = plot_energy_scan(
    #         scan_results["bond_lengths"],
    #         scan_results["energies"],
    #         title="O-H Bond Stretch in Water (HF/STO-3G)"
    #     )
    #     print(f"\nGenerated Base64 Image (first 50 chars): {img_b64[:50]}...")

    #     # THIS IS THE CRUCIAL PART TO TEST:
    #     html_file_path = "plot_demo.html"
    #     with open(html_file_path, "w") as f:
    #         f.write(f'<!DOCTYPE html><html><body><img src="data:image/png;base64,{img_b64}"/></body></html>')
    #     print(f"Image saved to {html_file_path}. Open this file in your web browser to view the plot.")

    # except Exception as e:
    #     print(f"Error during bond stretch/plot testing: {e}")


    # # (Original smiles_to_pyscf_input test can remain here too)
    # print("\n--- Original PySCF Input Test ---")
    # smiles = "CCO" # Ethanol
    # pyscf_input_code = smiles_to_xyz(smiles)
    # print("Generated PySCF xyz for Ethanol:")
    # print(pyscf_input_code)

    print("\n--- Testing 3D Molecular Visualization ---")
    try:
        # Example XYZ for water (can get from optimize_molecule result too)
        water_xyz = """5
Water molecule
C 0.000000 0.000000 0.000000
H -0.662308 0.852737 -0.084255
H -0.399859 -0.825892 -0.575263
H 0.083175 -0.291149 1.039820
H 0.978992 0.264303 -0.380303
"""
        html_viewer_content = generate_3d_viewer_html(water_xyz, "Water Molecule 3D View")
        
        html_file_path_3d = "molecule_3d_view.html"
        with open(html_file_path_3d, "w") as f:
            f.write(html_viewer_content)
        print(f"3D viewer HTML saved to {html_file_path_3d}. Opening in browser...")
        webbrowser.open(html_file_path_3d)

        # You could also chain it:
        optimized_methane = optimize_molecule(smiles_string="C")
        methane_3d_html = generate_3d_viewer_html(optimized_methane['optimized_geometry_xyz'], "Optimized Methane 3D View")
        with open("optimized_methane_3d.html", "w") as f:
            f.write(methane_3d_html)
        print(f"Optimized Methane 3D viewer HTML saved in /Users/lixinlu/Documents/mcp_server/optimized_methane_3d.html. Please open it in browser...")
        webbrowser.open("optimized_methane_3d.html")

    except Exception as e:
        print(f"Error during 3D visualization testing: {e}")


