"""
modules/pqr.py

Lightweight PQR utilities for the Brownian-dynamics scale.

Different tools emit slightly different PQR dialects: some include a chain-id
column, some append a trailing element column, and some (e.g. small-molecule
ligand PQRs written with one residue per atom) do neither. The reader here is
tolerant of all of these by anchoring the numeric fields from the right.

The writer slices a ``parmed.Structure`` down to a set of atoms (in topology
order), assigns MBondi2 radii, and emits a single PQR -- this is how a BD rigid
body's PQR is auto-generated from the molecular-dynamics system when the user
does not supply one.
"""

import os

from attrs import define

import parmed

import seekrflow.modules.base as base

try:
    import seekrtools.scripts.pqr_resid_for_each_atom as pqr_resid_for_each_atom
except ImportError:
    pqr_resid_for_each_atom = None


@define
class Pqr_atom:
    """
    A single parsed ATOM/HETATM record from a PQR file.
    """
    name: str
    resname: str
    resid: int
    x: float
    y: float
    z: float
    charge: float
    radius: float


def _is_float(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_pqr_atoms(filename: str) -> list[Pqr_atom]:
    """
    Parse the ATOM/HETATM records of a PQR file into an ordered list of
    ``Pqr_atom``.

    Tolerant of the two common PQR dialects: with or without a chain-id column,
    and with or without a trailing element column. The five numeric fields
    ``x y z charge radius`` are read from the right, and the residue id is the
    token immediately preceding them.
    """
    atoms: list[Pqr_atom] = []
    with open(filename, "r") as pqr_file:
        for line in pqr_file:
            tokens = line.split()
            if not tokens or tokens[0] not in ("ATOM", "HETATM"):
                continue
            # Drop a trailing element token if it is not numeric.
            if not _is_float(tokens[-1]):
                tokens = tokens[:-1]
            x, y, z, charge, radius = (float(t) for t in tokens[-5:])
            resid = int(tokens[-6])
            name = tokens[2]
            resname = tokens[3]
            atoms.append(Pqr_atom(
                name=name, resname=resname, resid=resid,
                x=x, y=y, z=z, charge=charge, radius=radius))
    return atoms


def write_pqr_from_parmed(
        parmed_complex: parmed.Structure,
        atom_indices: list[int],
        filename: str,
        one_resid_per_atom: bool = False,
        ) -> str:
    """
    Slice ``parmed_complex`` down to ``atom_indices`` (zero-based, in topology
    order), assign MBondi2 radii, and write a PQR to ``filename``.

    When ``one_resid_per_atom`` is True, rewrite the PQR so every atom gets its
    own residue id (required by some small-molecule BD setups). Returns the path
    that was written.
    """
    assert len(atom_indices) > 0, "Cannot write a PQR with no atoms."
    base.assign_mbondi2_radii_to_parmed_structure(parmed_complex)
    selection_string = "@" + ",".join(str(i + 1) for i in atom_indices)
    component_structure = parmed_complex[selection_string]
    component_structure.box = None
    if one_resid_per_atom:
        directory = os.path.dirname(filename)
        base_name = os.path.splitext(os.path.basename(filename))[0]
        tmp_pqr = os.path.join(directory, f"{base_name}_one_resid_tmp.pqr")
        component_structure.save(tmp_pqr, overwrite=True)
        assert pqr_resid_for_each_atom is not None, \
            "seekrtools is required to write per-atom-resid PQR files."
        pqr_resid_for_each_atom.pqr_resid_for_each_atom(tmp_pqr, filename)
    else:
        component_structure.save(filename, overwrite=True)
    return filename
