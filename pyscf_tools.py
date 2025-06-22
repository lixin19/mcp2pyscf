### python 3
# To generate sample pyscf input
#   can read in sample input from /sample/ directory
#   then add in geometry information
###

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, scf

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


# --- NEW FUNCTION FOR BOND STRETCHING ---
def run_bond_stretch_scan(smiles_string: str, atom1_idx: int, atom2_idx: int,
                          start_dist: float, end_dist: float, num_points: int) -> dict:
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

# Example usage for the new functions (for local testing)
if __name__ == "__main__":
    # Test bond stretching
    try:
        scan_results = run_bond_stretch_scan(
            smiles_string="O", # Water
            atom1_idx=0, atom2_idx=1, # O-H bond
            start_dist=0.8, end_dist=1.5, num_points=10
        )
        print("\nBond Scan Results:")
        print(scan_results)

        # Test plotting
        img_b64 = plot_energy_scan(
            scan_results["bond_lengths"],
            scan_results["energies"],
            title="O-H Bond Stretch in Water (HF/STO-3G)"
        )
        print(f"\nGenerated Base64 Image (first 50 chars): {img_b64[:50]}...")
        # To view locally, save this to an HTML file:
        # with open("plot_demo.html", "w") as f:
        #     f.write(f'<img src="data:image/png;base64,{img_b64}"/>')
        # print("Image saved to plot_demo.html. Open in browser to view.")

    except Exception as e:
        print(f"Error during bond stretch/plot testing: {e}")

    # (Original smiles_to_pyscf_input test can remain here too)
    print("\n--- Original PySCF Input Test ---")
    smiles = "CCO" # Ethanol
    pyscf_input_code = smiles_to_xyz(smiles)
    print("Generated PySCF xyz for Ethanol:")
    print(pyscf_input_code)




