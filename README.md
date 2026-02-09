# group-equivalent-ligand-atoms.py

## Description

This is an RDKit-based script that does the following:

1. Reads a molecular file (in MMCIF or PDB format) using Gemmi.
2. Finds ligand residues (non-water, non-standard amino acids, non-nucleotides).
3. Selects a set of uniquely named ligands - if there are more than one ligand with the same name, choose the ligand with the most atoms (or the first occuring if the sizes are tied).
3. For every ligand:
    * computes equivalence classes (that is, symmetry classes) for all the atoms using rdkit.Chem.CanonicalRankAtoms.
    * selects atoms from equivalence classes with more than one member
    * for each selected atom:
        * assigns same name to atoms from the same equivalence class, form the name from the chemical element and the equivalence class number separated by 'X'
        * prints a line with strings: {residue name} {original atom name} {assigned equivalence name}

## Setup to run

Create a Python virtual environment:

```bash
python3 -m venv myvenv
```

Activate a Python virtual environment:

```bash
source ./myvenv/bin/activate
```

Install RDKit and Gemmi:

```bash
pip install rdkit

pip install gemmi
```

Or, install using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

The script was tested with:

* Python version 3.12.3
* Gemmi version 0.7.4
* RDKIT version 2025.9.4

# Example

Get a molecular structure file, for example:

```bash
wget 'https://files.rcsb.org/download/2HHB.cif'
```

Then run the script and save the output to a file:

```bash
python3 ./group-equivalent-ligand-atoms.py ./2HHB.cif > ./table.tsv
```

The file `table.tsv` should contain a table where each row corresponds to an atom that belongs to a multi-member equivalence class within a ligand:

```
HEM    O1A    OX7
HEM    O2A    OX7
HEM    O1D    OX6
HEM    O2D    OX6
```

