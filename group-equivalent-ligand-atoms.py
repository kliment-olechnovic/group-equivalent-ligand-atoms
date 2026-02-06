#!/usr/bin/env python3

#
# DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
# Version 2, December 2004
#
# Copyright © 2026 Kliment Olechnovic
#
# Everyone is permitted to copy and distribute verbatim or modified
# copies of this license document, and changing it is allowed as long
# as the name is changed.
#
# DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
# TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
#
# 0. You just DO WHAT THE FUCK YOU WANT TO.
#

import sys
import os

import gemmi
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

AA = {
	"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
	"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}
NT = {"A","C","G","U","T","DA","DC","DG","DT","RA","RC","RG","RU"}

def is_polymer_resname(name: str) -> bool:
	name = (name or "").strip().upper()
	return name in AA or name in NT

def is_water(name: str) -> bool:
	return (name or "").strip().upper() in {"HOH", "WAT", "DOD"}

def normalize_ext(path: str) -> str:
	base = os.path.basename(path).lower()
	if base.endswith(".mmcif"):
		return ".cif"
	return os.path.splitext(base)[1]

def residue_to_pdb_block(res: gemmi.Residue, chain_id: str) -> str:
	st = gemmi.Structure()
	st.name = "one_res"
	st.cell = gemmi.UnitCell(1, 1, 1, 90, 90, 90)
	st.spacegroup_hm = "P 1"

	model = gemmi.Model("1")
	chain = gemmi.Chain(chain_id)

	r2 = gemmi.Residue()
	r2.name = res.name
	r2.seqid = gemmi.SeqId(res.seqid.num, res.seqid.icode)

	for a in res:
		a2 = gemmi.Atom()
		a2.name = a.name
		a2.element = a.element
		a2.pos = a.pos
		a2.occ = a.occ
		a2.b_iso = a.b_iso
		r2.add_atom(a2)

	chain.add_residue(r2)
	model.add_chain(chain)
	st.add_model(model)
	return st.make_pdb_string()

def rdkit_mol_from_residue(res: gemmi.Residue, chain_id: str) -> Chem.Mol:
	pdb_block = residue_to_pdb_block(res, chain_id)
	mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
	if mol is None:
		raise RuntimeError("RDKit failed to parse residue")

	try:
		rdDetermineBonds.DetermineBonds(mol)
	except Exception:
		rdDetermineBonds.DetermineConnectivity(mol)

	try:
		Chem.SanitizeMol(mol)
	except Exception:
		pass

	return mol

def equivalence_classes(mol: Chem.Mol):
	ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
	uniq = sorted(set(ranks))
	remap = {r:i for i, r in enumerate(uniq)}
	return [remap[r] for r in ranks]

def ligand_residues(structure: gemmi.Structure):
	for model in structure:
		for chain in model:
			chain_id = chain.name.strip() or "A"
			for res in chain:
				resname = (res.name or "").strip()
				if not resname:
					continue
				if is_water(resname):
					continue
				if is_polymer_resname(resname):
					continue
				if len(res) < 2:
					continue
				yield chain_id, res

def main():
	if len(sys.argv) != 2:
		print("Usage: python3 group-equivalent-ligand-atoms.py input.cif", file=sys.stderr)
		sys.exit(2)

	in_path = sys.argv[1]
	ext = normalize_ext(in_path)

	structure = gemmi.read_structure(in_path)

	print("residue_name\tatom_name\tequivalence_class")

	for chain_id, res in ligand_residues(structure):
		mol = rdkit_mol_from_residue(res, chain_id)
		eq = equivalence_classes(mol)

		for atom in mol.GetAtoms():
			info = atom.GetPDBResidueInfo()
			if info is None:
				continue
			resn = info.GetResidueName().strip()
			an   = info.GetName().strip()
			cls  = eq[atom.GetIdx()]
			print(f"{resn}\t{an}\t{cls}")

if __name__ == "__main__":
	main()
