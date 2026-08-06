"""
Membrane / lipid-atom helpers for equilibration restraints.

Ported from the legacy minimize_equilibrate script and aligned with
CHARMM-GUI Membrane Builder conventions (P Z-planar + lipid dihedrals).
"""

from __future__ import annotations

import math
import re
import typing

import numpy as np
import mdtraj

# Common lipid residue names in the CHARMM force field.
LIPID_NAMES = {
    "POPC", "POPE", "POPS", "POPG",
    "PSPC", "SDPC", "PAPC", "SSM",
    "DOPC", "DOPE", "DOPS", "DOPG",
    "DPPC", "DPPE", "DPPS", "DPPG",
    "DLPC", "DLPE", "DLPS", "DLPG",
    "DMPC", "DMPE", "DMPS", "DMPG",
    "CHOL", "CHL1",
    "PALM", "STEA", "OLEO", "LINO",
    "SAPI", "SAPI24", "SAPI25",
    "POPA", "DOPA", "DPPA", "DLPA", "DMPA",
}

CHOLESTEROL_NAMES = {"CHOL", "CHL1"}

# Target improper angle for glycerol C2 with atom order (C1, C3, C2, O21).
GLYCEROL_IMPROPER_THETA0_NEG = -2.0943951023931953  # -120 deg
GLYCEROL_IMPROPER_THETA0_POS = 2.0943951023931953   # +120 deg


def is_lipid_residue(residue: typing.Any) -> bool:
    """
    Whether a residue is a lipid (MEMB segment and/or known CHARMM resname).
    """
    if residue.is_protein or residue.is_water:
        return False
    if hasattr(residue, "segment_id") and residue.segment_id == "MEMB":
        return True
    if residue.name in LIPID_NAMES:
        return True
    return False


def _build_neighbor_map(traj: mdtraj.Trajectory) -> dict[int, dict[str, list[int]]]:
    neighbor_map: dict[int, dict[str, list[int]]] = {
        i: {"C": [], "H": [], "other": []} for i in range(traj.n_atoms)
    }
    for bond in traj.topology.bonds:
        atom1, atom2 = bond[0], bond[1]
        idx1, idx2 = atom1.index, atom2.index
        elem1, elem2 = atom1.element.symbol, atom2.element.symbol
        for a_idx, b_elem, b_idx in (
                (idx1, elem2, idx2), (idx2, elem1, idx1)):
            if b_elem == "C":
                neighbor_map[a_idx]["C"].append(b_idx)
            elif b_elem == "H":
                neighbor_map[a_idx]["H"].append(b_idx)
            else:
                neighbor_map[a_idx]["other"].append(b_idx)
    return neighbor_map


def _dihedral_radians(
        positions: np.ndarray,
        indices: tuple[int, int, int, int],
        ) -> float:
    """Compute a torsion angle in radians from Cartesian positions (nm)."""
    p0, p1, p2, p3 = (positions[i] for i in indices)
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm
    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)
    return float(math.atan2(y, x))


def glycerol_improper_theta0(
        positions: np.ndarray,
        improper: tuple[int, int, int, int],
        ) -> float:
    """
    Choose ±120° for a glycerol improper based on the structure.

    Atom order is (C1, C3, C2, O21); standard sn lipids prefer −120°.
    """
    theta = _dihedral_radians(positions, improper)
    if abs(theta - GLYCEROL_IMPROPER_THETA0_POS) < abs(
            theta - GLYCEROL_IMPROPER_THETA0_NEG):
        return GLYCEROL_IMPROPER_THETA0_POS
    return GLYCEROL_IMPROPER_THETA0_NEG


def identify_lipid_atoms(
        solvated_mdtraj: mdtraj.Trajectory,
        ) -> dict[str, typing.Any]:
    """
    Identify lipid atom groups for membrane equilibration restraints.

    Returns keys:
      * ``phosphorus`` — P atom indices (lipid residues only when possible)
      * ``cis_double_bonds`` — list of 4-tuples (C–C=C–C)
      * ``glycerol_impropers`` — list of 4-tuples (C1, C3, C2, O21)
    """
    result: dict[str, typing.Any] = {
        "phosphorus": [],
        "cis_double_bonds": [],
        "glycerol_impropers": [],
    }

    print("Detecting lipid atoms using CHARMM naming conventions...")
    try:
        p_all = list(solvated_mdtraj.topology.select("element P"))
    except Exception as exc:
        print(f"WARNING: Phosphorus detection failed: {exc}")
        print("Proceeding without lipid restraints.")
        return result

    p_indices = [
        int(i) for i in p_all
        if is_lipid_residue(solvated_mdtraj.topology.atom(int(i)).residue)
    ]
    # Fall back to all P if lipid filtering removes everything (odd topologies).
    if not p_indices and p_all:
        p_indices = [int(i) for i in p_all]
    if not p_indices:
        print("WARNING: No phosphorus atoms found. Skipping lipid restraints.")
        return result

    print(f"Found {len(p_indices)} phosphorus atoms in system.")
    result["phosphorus"] = p_indices

    neighbor_map = _build_neighbor_map(solvated_mdtraj)

    tail_pattern = re.compile(r"^C[23]\d+$")
    tail_carbons = [
        atom.index for atom in solvated_mdtraj.topology.atoms
        if is_lipid_residue(atom.residue)
        and tail_pattern.match(atom.name)
        and atom.element.symbol == "C"
    ]
    print(f"Found {len(tail_carbons)} tail carbon atoms (C2X/C3X pattern).")

    sp2_carbons = []
    for carbon_idx in tail_carbons:
        neighbors = neighbor_map[carbon_idx]
        if len(neighbors["C"]) == 2 and len(neighbors["H"]) == 1:
            sp2_carbons.append(carbon_idx)
    print(f"Identified {len(sp2_carbons)} sp2 carbons (potential double bonds).")

    double_bond_pairs: list[tuple[int, int]] = []
    sp2_set = set(sp2_carbons)
    for sp2_idx in sp2_carbons:
        for neighbor_idx in neighbor_map[sp2_idx]["C"]:
            if neighbor_idx in sp2_set:
                pair = tuple(sorted((sp2_idx, neighbor_idx)))
                if pair not in double_bond_pairs:
                    double_bond_pairs.append(pair)  # type: ignore[arg-type]
    print(f"Found {len(double_bond_pairs)} C=C double bonds.")

    for c_i, c_j in double_bond_pairs:
        c_i_neighbors = [n for n in neighbor_map[c_i]["C"] if n != c_j]
        c_j_neighbors = [n for n in neighbor_map[c_j]["C"] if n != c_i]
        if c_i_neighbors and c_j_neighbors:
            result["cis_double_bonds"].append(
                (c_i_neighbors[0], c_i, c_j, c_j_neighbors[0]))
    print(
        f"Created {len(result['cis_double_bonds'])} cis double-bond "
        "dihedral restraints.")

    residue_atoms: dict[tuple[int, int], dict] = {}
    for atom in solvated_mdtraj.topology.atoms:
        if (is_lipid_residue(atom.residue)
                and atom.residue.name not in CHOLESTEROL_NAMES):
            res_key = (atom.residue.chain.index, atom.residue.index)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = {
                    "atoms": [], "residue": atom.residue}
            residue_atoms[res_key]["atoms"].append(atom)

    print(f"Processing {len(residue_atoms)} non-cholesterol lipid residues...")
    impropers_found = 0
    for res_data in residue_atoms.values():
        atoms = res_data["atoms"]
        by_name = {atom.name: atom for atom in atoms}
        c1_atom = by_name.get("C1")
        c2_atom = by_name.get("C2")
        c3_atom = by_name.get("C3")
        o21_atom = by_name.get("O21")
        if c1_atom is None or c2_atom is None or c3_atom is None \
                or o21_atom is None:
            continue
        c2_neighbors = (
            neighbor_map[c2_atom.index]["C"]
            + neighbor_map[c2_atom.index]["other"])
        if not (c1_atom.index in c2_neighbors
                and c3_atom.index in c2_neighbors
                and o21_atom.index in c2_neighbors):
            continue
        result["glycerol_impropers"].append(
            (c1_atom.index, c3_atom.index, c2_atom.index, o21_atom.index))
        impropers_found += 1

    print(
        f"Created {impropers_found} glycerol improper dihedral restraints "
        "(C1-C3-C2-O21).")
    return result
