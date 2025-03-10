"""
This script scans a directory for all PDB files and compares them to two reference ADK structures (1AKE - closed, 4AKE - open).
It produces per-structure RMSDs, then visualizes the results in a scatter plot, X-axis = RMSD to 1AKE, Y-axis = RMSD to 4AKE.

How to use:
    python pdb_compare.py --pdb_folder /path/to/folder
"""

import os
import re
import argparse
import requests
import matplotlib.pyplot as plt
from datetime import datetime

from adjustText import adjust_text
from Bio import PDB
from Bio.PDB import PDBParser, Superimposer

# HELPER FUNCTIONS

def fetch_pdb(pdb_id, save_dir="."):
    #:param pdb_id (e.g. '1AKE')
    #:param save_dir: Directory to save the downloaded file
    #:return: Path to the fetched PDB file
    
    pdb_path = os.path.join(save_dir, f"{pdb_id}.pdb")
    if not os.path.isfile(pdb_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Downloading {pdb_id}.pdb from RCSB...")
        r = requests.get(url)
        if r.status_code == 200:
            with open(pdb_path, 'w') as f:
                f.write(r.text)
        else:
            raise IOError(f"Could not fetch PDB {pdb_id} from RCSB (status code={r.status_code}).")
    return pdb_path


def get_chain_residues(structure, chain_id='A'):
    #:param structure: biopython structure object
    #:param chain_id: the chain you want to extract
    #:return: list of residues in the specified chain
    
    model = structure[0]  # only dealing with first model
    chain = model[chain_id] # only dealing with chain A since AF2 gives only chain A
    residues = [res for res in chain]
    return residues


def get_atoms_for_alignment(residues, atom_type="CA"):
    # Given a list of residues, extract specific atoms (e.g., alpha carbons) for alignment.
    #:param residues: list of Biopython Residue objects
    #:param atom_type: which atom(s) to extract, default CA for alpha-carbons
    #:return: list of Biopython Atom objects
    
    atoms = []
    for residue in residues:
        if atom_type in residue:
            atoms.append(residue[atom_type])
    return atoms


def compute_rmsd(moving_struct, ref_struct, chain_id='A'):
    # Compute RMSD by aligning chain 'A' of the moving_struct onto
    # chain 'A' of the ref_struct using alpha carbons.
    # :param moving_struct: Biopython structure object to align (the "moving" target)
    # :param ref_struct: Biopython structure object that is the reference
    # :param chain_id: the chain ID to align (assumes same chain ID in both)
    # :return: RMSD (float)
    
    # Extract chain residues
    moving_residues = get_chain_residues(moving_struct, chain_id=chain_id)
    ref_residues = get_chain_residues(ref_struct, chain_id=chain_id)

    # Extract CA atoms
    moving_atoms = get_atoms_for_alignment(moving_residues, atom_type="CA")
    ref_atoms = get_atoms_for_alignment(ref_residues, atom_type="CA")

    # If not the same length, a simple approach is to align up to the min length
    if len(moving_atoms) != len(ref_atoms):
        min_len = min(len(moving_atoms), len(ref_atoms))
        moving_atoms = moving_atoms[:min_len]
        ref_atoms = ref_atoms[:min_len]

    # Superimpose
    sup = Superimposer()
    sup.set_atoms(ref_atoms, moving_atoms) # set atoms to align
    sup.apply(moving_struct.get_atoms())  # apply the transformation to the moving structure
    return sup.rms


# -------------------------------------------------------------
# MAIN LOGIC
# -------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compare ADK PDB structures to 1AKE & 4AKE")
    parser.add_argument("--pdb_folder", required=True,
                        help="Path to the folder containing PDB files to compare.")
    parser.add_argument("--chain_id", default="A",
                        help="Chain ID to align on. Default='A'")
    args = parser.parse_args()

    pdb_folder = args.pdb_folder
    chain_id = args.chain_id

    # 1. Gather local PDB files from the given folder
    print(f"Scanning {pdb_folder} for .pdb files...")
    all_pdb_files = []
    for root, dirs, files in os.walk(pdb_folder):
        for f in files:
            if f.endswith(".pdb"):
                all_pdb_files.append(os.path.join(root, f))

    if not all_pdb_files:
        print("No PDB files found in the specified folder.")
        return

    # 2. Download/reference the known ADK structures: 1AKE (closed), 4AKE (open)
    ref_ids = ["1AKE", "4AKE"]
    references = {}
    for r_id in ref_ids:
        references[r_id] = fetch_pdb(r_id, save_dir=pdb_folder)

    # 3. Parse reference structures
    parser_biopdb = PDBParser(QUIET=True)
    ref_structs = {}
    for r_id in ref_ids:
        ref_structs[r_id] = parser_biopdb.get_structure(r_id, references[r_id])

    # 4. Compute RMSD for each local file vs 1AKE and 4AKE
    results = []
    for pdb_file in all_pdb_files:
        structure_id = os.path.basename(pdb_file)
        moving_structure = parser_biopdb.get_structure(structure_id, pdb_file)

        # Align to 1AKE (closed), use .copy to preserve the reference structure
        rmsd_1AKE = compute_rmsd(moving_structure.copy(), ref_structs["1AKE"].copy(), chain_id=chain_id)

        # Align to 4AKE (open)
        rmsd_4AKE = compute_rmsd(moving_structure.copy(), ref_structs["4AKE"].copy(), chain_id=chain_id)

        results.append((structure_id, rmsd_1AKE, rmsd_4AKE))

    # 5. Print results
    print("\nComparison Results (RMSD in Å):")
    max_name_length = max(len(s_id) for s_id, _, _ in results) if results else 40
    column_width = max(40, max_name_length + 4)  # Ensure at least 40 chars

    # Print header and row with dynamic width
    print(f"|{'Structure':{column_width}s}|{'RMSD vs 1AKE':>12s}|{'RMSD vs 4AKE':>12s}|")
    for (s_id, rms1, rms4) in results:
        print(f"|{s_id:{column_width}s}|{rms1:12.3f}|{rms4:12.3f}|")

    # 6. SCATTER PLOT and save to results:
    #    x-axis = RMSD vs 1AKE
    #    y-axis = RMSD vs 4AKE
    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)  

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file_name = os.path.basename(pdb_folder)
    output_file = os.path.join(results_dir, f"adk_rmsd_scatter_{file_name}_{timestamp}.png")

    # Extract structure names and RMSD values
    structure_labels = [r[0] for r in results]
    rmsd_1 = [r[1] for r in results]  # RMSD vs 1AKE
    rmsd_4 = [r[2] for r in results]  # RMSD vs 4AKE

    plt.figure(figsize=(8, 6))

    # Scatter plot with improved marker properties
    plt.scatter(rmsd_1, rmsd_4, c='blue', edgecolors='black', alpha=0.7, linewidths=0.5)

    plt.xlabel("RMSD vs 1AKE (closed)")
    plt.ylabel("RMSD vs 4AKE (open)")
    plt.title("Distribution of AF2 ADK Predictions: RMSD to 1AKE vs 4AKE")

    # # Optionally label points with the structure ID
    # texts = []
    # for i, label in enumerate(structure_labels):
    #     texts.append(plt.text(rmsd_1[i], rmsd_4[i], label, fontsize=8, ha='left', va='bottom'))

    # # Adjust text to prevent overlap
    # adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Scatter plot saved as {output_file}")
    plt.show()

if __name__ == "__main__":
    main()
