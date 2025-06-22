### python 3
# To generate geometry from SMILES using rdkit
###

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


