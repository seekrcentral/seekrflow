"""
modules/pqr.py

Lightweight PQR utilities for the Brownian-dynamics scale.

Different tools emit slightly different PQR dialects: some include a chain-id
column, some append a trailing element column, and some (e.g. small-molecule
ligand PQRs written with one residue per atom) do neither. The reader here is
tolerant of all of these by anchoring the numeric fields from the right.

The writer emits a strict fixed-column PQR (ATOM records only, no chain id,
resnames truncated to 3 characters) from either a list of ``Pqr_atom`` or a
sliced ``parmed.Structure``. This is how a BD rigid body's PQR is auto-generated
from the molecular-dynamics system when the user does not supply one.
"""

from __future__ import annotations

import os

from attrs import define

import parmed

import seekrflow.modules.base as base

# Strict classic PQR line without chain / element columns:
# ATOM  serial name resname resid    x       y       z      charge  radius
_PQR_ATOM_FMT = (
    "ATOM  {serial:5d} {name:4s} {resname:3s} {resid:5d}    "
    "{x:8.3f}{y:8.3f}{z:8.3f} {charge:7.4f} {radius:7.4f}\n"
)


@define
class Pqr_atom:
    """
    A single ATOM record for PQR read/write.
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


def _format_atom_name(name: str) -> str:
    """Pad/truncate atom name to the 4-column PDB/PQR name field."""
    text = (name or "").strip()
    if len(text) > 4:
        text = text[:4]
    return f"{text:<4s}"


def _format_resname(resname: str) -> str:
    """Truncate residue name to 3 characters (left) for classic PQR columns."""
    text = (resname or "").strip()
    if len(text) > 3:
        text = text[:3]
    return f"{text:>3s}"


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


def write_pqr_atoms(atoms: list[Pqr_atom], filename: str) -> str:
    """
    Write ``atoms`` as a strict fixed-column PQR (ATOM only, no chain id).

    Residue names longer than 3 characters are truncated on the left
    (``CHL1`` → ``CHL``). Returns ``filename``.
    """
    if not atoms:
        raise ValueError("Cannot write a PQR with no atoms.")
    directory = os.path.dirname(filename)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(filename, "w") as handle:
        for index, atom in enumerate(atoms):
            serial = index + 1
            if serial > 99999:
                # Classic 5-column serial field overflows; keep writing with
                # wider serial so large BD bodies are not silently corrupted.
                line = (
                    f"ATOM  {serial:d} {_format_atom_name(atom.name)} "
                    f"{_format_resname(atom.resname)} {int(atom.resid):5d}    "
                    f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f} "
                    f"{atom.charge:7.4f} {atom.radius:7.4f}\n"
                )
            else:
                line = _PQR_ATOM_FMT.format(
                    serial=serial,
                    name=_format_atom_name(atom.name),
                    resname=_format_resname(atom.resname),
                    resid=int(atom.resid),
                    x=float(atom.x),
                    y=float(atom.y),
                    z=float(atom.z),
                    charge=float(atom.charge),
                    radius=float(atom.radius),
                )
            handle.write(line)
    return filename


def structure_to_pqr_atoms(
        parmed_structure: parmed.Structure,
        atom_indices: list[int],
        *,
        one_resid_per_atom: bool = False,
        ) -> list[Pqr_atom]:
    """
    Build ``Pqr_atom`` records from ``parmed_structure`` for ``atom_indices``
    (zero-based, topology order).

    Uses each atom's coordinates, charge, and ``solvent_radius`` (assign MBondi2
    before calling). Does not use ParmEd selection or ``Structure.save``.
    """
    n_atoms = len(parmed_structure.atoms)
    atoms: list[Pqr_atom] = []
    for out_index, atom_index in enumerate(atom_indices):
        if atom_index < 0 or atom_index >= n_atoms:
            raise IndexError(
                f"atom index {atom_index} out of range for structure with "
                f"{n_atoms} atoms")
        atom = parmed_structure.atoms[atom_index]
        resid = (out_index + 1) if one_resid_per_atom else int(atom.residue.number)
        radius = getattr(atom, "solvent_radius", None)
        if radius is None:
            radius = getattr(atom, "rmin", 0.0) or 0.0
        atoms.append(Pqr_atom(
            name=str(atom.name),
            resname=str(atom.residue.name),
            resid=int(resid),
            x=float(atom.xx),
            y=float(atom.xy),
            z=float(atom.xz),
            charge=float(atom.charge),
            radius=float(radius),
        ))
    return atoms


def write_pqr_from_parmed(
        parmed_complex: parmed.Structure,
        atom_indices: list[int],
        filename: str,
        one_resid_per_atom: bool = False,
        ) -> str:
    """
    Slice ``parmed_complex`` down to ``atom_indices`` (zero-based, in topology
    order), assign MBondi2 radii, and write a strict classic PQR to ``filename``.

    When ``one_resid_per_atom`` is True, each atom gets residue id 1..N (required
    by some small-molecule BD setups). Returns the path that was written.
    """
    assert len(atom_indices) > 0, "Cannot write a PQR with no atoms."
    base.assign_mbondi2_radii_to_parmed_structure(parmed_complex)
    atoms = structure_to_pqr_atoms(
        parmed_complex,
        atom_indices,
        one_resid_per_atom=one_resid_per_atom,
    )
    return write_pqr_atoms(atoms, filename)
