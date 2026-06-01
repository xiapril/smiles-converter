#!/usr/bin/env python3
"""
SMILES -> Oligo Notation + Molecular Weight + Conjugate Detection

Outputs from ONE SMILES input (single strand or duplex with "."):
  1) Thermo BioPharma Finder Sequence (FULL style)
  2) Agilent BioConfirm Sequence (unchanged)
  3) Molecular Weight Summary (full product, each strand)
  4) Conjugate Detection with formula + MW

Usage:
  python smiles_to_oligo_notation.py "YOUR_SMILES_HERE"
  python smiles_to_oligo_notation.py        (will prompt for input)

Requires: RDKit (rdkit)
"""

import sys
from collections import Counter
from rdkit import Chem
from rdkit.Chem import Descriptors


# ============================================================
# Atomic masses for fragment MW calculation
# ============================================================
MONO_MASS = {
    "H": 1.00783, "C": 12.00000, "N": 14.00307, "O": 15.99491,
    "F": 18.99840, "P": 30.97376, "S": 31.97207, "Cl": 34.96885,
    "Br": 78.91834, "I": 126.90447
}
AVG_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "P": 30.974, "S": 32.065, "Cl": 35.453,
    "Br": 79.904, "I": 126.904
}


# ============================================================
# Helper: methyl / ethyl detection
# ============================================================
def _is_methyl_carbon(mol, c_idx, attached_to_idx=None):
    a = mol.GetAtomWithIdx(c_idx)
    if a.GetSymbol() != "C":
        return False
    heavy = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() != "H"]
    if attached_to_idx is None:
        return len(heavy) == 1
    return len(heavy) == 1 and heavy[0] == attached_to_idx


def _is_ethyl_from_oxygen(mol, o_idx, p_idx, exclude=set()):
    o = mol.GetAtomWithIdx(o_idx)
    others = [n for n in o.GetNeighbors() if n.GetIdx() != p_idx]
    if len(others) != 1:
        return False
    c1 = others[0]
    if c1.GetSymbol() != "C" or c1.GetIdx() in exclude:
        return False
    c1_idx = c1.GetIdx()
    c1_nbrs = [n for n in c1.GetNeighbors() if n.GetSymbol() != "H" and n.GetIdx() != o_idx]
    if len(c1_nbrs) != 1:
        return False
    c2 = c1_nbrs[0]
    if c2.GetSymbol() != "C" or c2.GetIdx() in exclude:
        return False
    return _is_methyl_carbon(mol, c2.GetIdx(), attached_to_idx=c1_idx)


# ============================================================
# Fragment formula + MW computation
# ============================================================
def compute_fragment_formula(mol, atom_indices):
    atom_set = set(atom_indices)
    counts = Counter()
    for idx in atom_set:
        atom = mol.GetAtomWithIdx(idx)
        counts[atom.GetSymbol()] += 1
        counts["H"] += atom.GetNumImplicitHs()
    for idx in atom_set:
        for bond in mol.GetAtomWithIdx(idx).GetBonds():
            if bond.GetOtherAtomIdx(idx) not in atom_set:
                counts["H"] += 1
    formula = ""
    if "C" in counts:
        formula += f"C{counts['C']}" if counts['C'] > 1 else "C"
        del counts["C"]
    if "H" in counts:
        formula += f"H{counts['H']}" if counts['H'] > 1 else "H"
        del counts["H"]
    for sym in sorted(counts.keys()):
        formula += f"{sym}{counts[sym]}" if counts[sym] > 1 else sym
    return formula


def compute_fragment_mw(mol, atom_indices):
    atom_set = set(atom_indices)
    counts = Counter()
    for idx in atom_set:
        atom = mol.GetAtomWithIdx(idx)
        counts[atom.GetSymbol()] += 1
        counts["H"] += atom.GetNumImplicitHs()
    for idx in atom_set:
        for bond in mol.GetAtomWithIdx(idx).GetBonds():
            if bond.GetOtherAtomIdx(idx) not in atom_set:
                counts["H"] += 1
    mono = sum(MONO_MASS.get(s, 0) * c for s, c in counts.items())
    avg = sum(AVG_MASS.get(s, 0) * c for s, c in counts.items())
    return mono, avg


def _find_phosphorus_via_oxygen(mol, carbon_idx, ring_set):
    if carbon_idx is None:
        return None
    for nbr in mol.GetAtomWithIdx(carbon_idx).GetNeighbors():
        if nbr.GetSymbol() == "O" and nbr.GetIdx() not in ring_set:
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetSymbol() == "P":
                    return nbr2.GetIdx()
    return None

# ============================================================
# Base extraction + Thermo base symbols
# ============================================================
def _find_base_atoms(mol, gly_n, ring_set):
    base_atoms = set()
    stack = [gly_n]
    while stack:
        curr = stack.pop()
        if curr in base_atoms or curr in ring_set:
            continue
        base_atoms.add(curr)
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx in base_atoms or nidx in ring_set or nbr.GetSymbol() == "P":
                continue
            stack.append(nidx)
    return base_atoms


def _pick_pyrimidine_ring(mol, base_atoms):
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 6 and all(r in base_atoms for r in ring):
            if sum(1 for r in ring if mol.GetAtomWithIdx(r).GetSymbol() == "N") >= 2:
                return set(ring)
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 6 and any(r in base_atoms for r in ring):
            return set(ring)
    return set()


def _pyrimidine_has_methyl(mol, ring_set, base_atoms):
    for ridx in ring_set:
        a = mol.GetAtomWithIdx(ridx)
        if a.GetSymbol() != "C":
            continue
        for nbr in a.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx in ring_set:
                continue
            if nidx in base_atoms and nbr.GetSymbol() == "C":
                if _is_methyl_carbon(mol, nidx, attached_to_idx=ridx):
                    return True
    return False


def _pyrimidine_has_methoxy(mol, ring_set, base_atoms):
    for ridx in ring_set:
        a = mol.GetAtomWithIdx(ridx)
        if a.GetSymbol() != "C":
            continue
        for nbr in a.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx in ring_set:
                continue
            if nidx in base_atoms and nbr.GetSymbol() == "O":
                o = mol.GetAtomWithIdx(nidx)
                other = [x for x in o.GetNeighbors() if x.GetIdx() != ridx]
                if len(other) == 1 and other[0].GetSymbol() == "C":
                    if _is_methyl_carbon(mol, other[0].GetIdx(), attached_to_idx=nidx):
                        return True
    return False


def _pyrimidine_has_thio(mol, base_atoms):
    return any(mol.GetAtomWithIdx(i).GetSymbol() == "S" for i in base_atoms)


def _pyrimidine_has_ethyl_like(mol, ring_set, base_atoms):
    for ridx in ring_set:
        a = mol.GetAtomWithIdx(ridx)
        if a.GetSymbol() != "C":
            continue
        for nbr in a.GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx in ring_set:
                continue
            if nidx in base_atoms and nbr.GetSymbol() == "C":
                c1_nbrs = [n for n in mol.GetAtomWithIdx(nidx).GetNeighbors()
                           if n.GetSymbol() != "H" and n.GetIdx() != ridx]
                if len(c1_nbrs) == 1 and c1_nbrs[0].GetSymbol() == "C":
                    if _is_methyl_carbon(mol, c1_nbrs[0].GetIdx(), attached_to_idx=nidx):
                        return True
    return False


def _purine_has_exocyclic_amine(mol, base_atoms):
    for i in base_atoms:
        a = mol.GetAtomWithIdx(i)
        if a.GetSymbol() != "N" or a.GetIsAromatic():
            continue
        heavy = [n for n in a.GetNeighbors() if n.GetSymbol() != "H"]
        if len(heavy) == 1 and heavy[0].GetIdx() in base_atoms and heavy[0].GetSymbol() == "C":
            return True
    return False


def _purine_has_methyl_on_ring(mol, base_atoms):
    for i in base_atoms:
        a = mol.GetAtomWithIdx(i)
        if a.GetSymbol() != "N" or not a.GetIsAromatic():
            continue
        for nbr in a.GetNeighbors():
            if nbr.GetSymbol() == "C" and nbr.GetIdx() in base_atoms:
                if _is_methyl_carbon(mol, nbr.GetIdx(), attached_to_idx=i):
                    return True
    return False


def identify_thermo_base(mol, gly_n, ring_set):
    base_atoms = _find_base_atoms(mol, gly_n, ring_set)
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in base_atoms]
    o_count = syms.count("O")
    is_purine = len(base_atoms) >= 10

    if not is_purine:
        ring6 = _pick_pyrimidine_ring(mol, base_atoms)
        has_methyl = _pyrimidine_has_methyl(mol, ring6, base_atoms) if ring6 else False
        has_methoxy = _pyrimidine_has_methoxy(mol, ring6, base_atoms) if ring6 else False
        has_thio = _pyrimidine_has_thio(mol, base_atoms)
        has_eth = _pyrimidine_has_ethyl_like(mol, ring6, base_atoms) if ring6 else False

        if o_count >= 2:
            if has_thio:
                return "B"
            if has_methoxy:
                return "D"
            if has_eth:
                return "Q"
            if has_methyl:
                return "T"
            return "U"
        if o_count == 1:
            return "S" if has_methyl else "C"
        return "?"

    if o_count == 0:
        return "A"
    if _purine_has_methyl_on_ring(mol, base_atoms):
        return "N"
    if _purine_has_exocyclic_amine(mol, base_atoms):
        return "G"
    return "I"


# ============================================================
# Thermo 2' ribose symbol
# ============================================================
def identify_thermo_2prime(mol, c2_idx, ring_set):
    if c2_idx is None:
        return "?"
    exo = [(n.GetSymbol(), n.GetIdx()) for n in mol.GetAtomWithIdx(c2_idx).GetNeighbors()
           if n.GetIdx() not in ring_set]
    if any(sym == "F" for sym, _ in exo):
        return "f"
    o_idx = None
    for sym, idx in exo:
        if sym == "O":
            o_idx = idx
            break
    if o_idx is None:
        return "d"
    o = mol.GetAtomWithIdx(o_idx)
    o_others = [n for n in o.GetNeighbors() if n.GetIdx() != c2_idx]
    if len(o_others) == 0:
        return "r"
    if len(o_others) == 1 and o_others[0].GetSymbol() == "C":
        c1_idx = o_others[0].GetIdx()
        if _is_methyl_carbon(mol, c1_idx, attached_to_idx=o_idx):
            return "m"
        visited = {c2_idx, o_idx}
        queue = [c1_idx]
        heavy = []
        hetero = {"O": 0, "N": 0, "S": 0}
        while queue:
            cur = queue.pop(0)
            if cur in visited:
                continue
            visited.add(cur)
            a = mol.GetAtomWithIdx(cur)
            sym = a.GetSymbol()
            if sym != "H":
                heavy.append(cur)
                if sym in hetero:
                    hetero[sym] += 1
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetSymbol() != "P":
                    queue.append(nbr.GetIdx())
        chain_len = len(heavy)
        if hetero["S"] >= 1 and chain_len >= 6:
            return "a"
        if hetero["N"] >= 1 and chain_len >= 3:
            return "n"
        # MOE check
        for nbr_c2 in mol.GetAtomWithIdx(c1_idx).GetNeighbors():
            if nbr_c2.GetSymbol() != "C" or nbr_c2.GetIdx() == o_idx:
                continue
            c2chain = nbr_c2.GetIdx()
            for nbr_o2 in mol.GetAtomWithIdx(c2chain).GetNeighbors():
                if nbr_o2.GetSymbol() != "O" or nbr_o2.GetIdx() == o_idx:
                    continue
                o2 = mol.GetAtomWithIdx(nbr_o2.GetIdx())
                o2_others = [x for x in o2.GetNeighbors() if x.GetIdx() != c2chain]
                if len(o2_others) == 1 and o2_others[0].GetSymbol() == "C":
                    if _is_methyl_carbon(mol, o2_others[0].GetIdx(), attached_to_idx=nbr_o2.GetIdx()):
                        return "e"
        if chain_len >= 12 and (hetero["O"] + hetero["N"] + hetero["S"]) <= 1:
            return "h"
        return "m"
    return "?"


# ============================================================
# Backbone linker classification (Thermo)
# ============================================================
def classify_thermo_linker(mol, p_idx, exclude=set()):
    if p_idx is None:
        return None
    p = mol.GetAtomWithIdx(p_idx)
    nbr_syms = [n.GetSymbol() for n in p.GetNeighbors()]
    if "P" in nbr_syms:
        return "t"
    if "N" in nbr_syms:
        return "n"
    has_s = "S" in nbr_syms
    has_ethyl = False
    for nbr in p.GetNeighbors():
        if nbr.GetSymbol() == "O":
            if _is_ethyl_from_oxygen(mol, nbr.GetIdx(), p_idx, exclude=exclude):
                has_ethyl = True
                break
    if has_ethyl and has_s:
        return "x"
    if has_ethyl:
        return "y"
    if has_s:
        return "s"
    return "p"


def is_phosphorothioate(mol, p_idx):
    if p_idx is None:
        return False
    return any(n.GetSymbol() == "S" for n in mol.GetAtomWithIdx(p_idx).GetNeighbors())

# ============================================================
# 5' vinyl phosphate detection
# ============================================================
def detect_5prime_vp(mol, nuc):
    c5_idx = nuc.get("c5")
    if c5_idx is None:
        return False
    c5 = mol.GetAtomWithIdx(c5_idx)
    for nbr in c5.GetNeighbors():
        if nbr.GetSymbol() != "C" or nbr.GetIdx() == nuc.get("c4", -1):
            continue
        if any(n2.GetSymbol() == "P" for n2 in nbr.GetNeighbors()):
            return True
        for n2 in nbr.GetNeighbors():
            if n2.GetSymbol() == "C":
                if any(n3.GetSymbol() == "P" for n3 in n2.GetNeighbors()):
                    return True
    return False


# ============================================================
# Conjugate detection
# ============================================================
def detect_2prime_conjugate(mol, c2_idx, ring_set):
    if c2_idx is None:
        return None
    o_idx = None
    for nbr in mol.GetAtomWithIdx(c2_idx).GetNeighbors():
        if nbr.GetIdx() not in ring_set and nbr.GetSymbol() == "O":
            o_idx = nbr.GetIdx()
            break
    if o_idx is None:
        return None
    o_others = [n for n in mol.GetAtomWithIdx(o_idx).GetNeighbors() if n.GetIdx() != c2_idx]
    if len(o_others) == 0:
        return None
    if len(o_others) == 1 and o_others[0].GetSymbol() == "C":
        if _is_methyl_carbon(mol, o_others[0].GetIdx(), attached_to_idx=o_idx):
            return None
    visited = {c2_idx}
    queue = [o_idx]
    conj = set()
    hetero = {"O": 0, "N": 0, "S": 0}
    while queue:
        cur = queue.pop(0)
        if cur in visited:
            continue
        visited.add(cur)
        a = mol.GetAtomWithIdx(cur)
        if a.GetSymbol() == "P":
            continue
        conj.add(cur)
        if a.GetSymbol() in hetero:
            hetero[a.GetSymbol()] += 1
        for nb in a.GetNeighbors():
            if nb.GetIdx() not in visited:
                queue.append(nb.GetIdx())
    if len(conj) <= 5:
        return None
    desc = "unknown conjugate"
    if hetero["S"] >= 1:
        desc = "ADS-like (sulfur-containing)"
    elif len(conj) >= 12 and (hetero["O"] + hetero["N"] + hetero["S"]) <= 1:
        desc = "HD-like (long alkyl)"
    elif hetero["N"] >= 1:
        desc = "amino-linked conjugate"
    return conj, desc


def detect_3prime_conjugate(mol, p3_idx, all_nuc_atoms, all_p5):
    if p3_idx is None or p3_idx in all_p5:
        return None
    visited = set(all_nuc_atoms)
    visited.discard(p3_idx)
    queue = [p3_idx]
    conj = set()
    while queue:
        cur = queue.pop(0)
        if cur in visited:
            continue
        visited.add(cur)
        conj.add(cur)
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            if nb.GetIdx() not in visited:
                queue.append(nb.GetIdx())
    non_po = [i for i in conj if mol.GetAtomWithIdx(i).GetSymbol() not in ("P", "O", "H")]
    if len(non_po) < 6:
        return None
    desc = "3' terminal conjugate"
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in conj]
    if syms.count("N") >= 3 and len(conj) > 30:
        desc = "3' terminal conjugate (GalNAc-like)"
    return conj, desc


def detect_5prime_conjugate(mol, nuc, all_nuc_atoms):
    p5 = nuc.get("p5")
    if p5 is None:
        return None
    visited = set(all_nuc_atoms)
    visited.discard(p5)
    queue = [p5]
    conj = set()
    while queue:
        cur = queue.pop(0)
        if cur in visited:
            continue
        visited.add(cur)
        conj.add(cur)
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            if nb.GetIdx() not in visited:
                queue.append(nb.GetIdx())
    non_po = [i for i in conj if mol.GetAtomWithIdx(i).GetSymbol() not in ("P", "O", "H")]
    if len(non_po) < 6:
        return None
    return conj, "5' terminal conjugate"


# ============================================================
# Parse one strand -> nucleotides, order, VP, conjugates, MW
# ============================================================
def parse_strand(smiles_str):
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        raise ValueError("RDKit could not parse the SMILES string.")

    mono_mw = Descriptors.ExactMolWt(mol)
    avg_mw = Descriptors.MolWt(mol)

    ring_info = mol.GetRingInfo()
    nucleotides = []
    all_nuc_atoms = set()
    sugar_like = set()

    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if syms.count("O") != 1 or syms.count("C") != 4:
            continue
        ring_set = set(ring)
        sugar_like |= ring_set
        ring_o = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
        ring_o_nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(ring_o).GetNeighbors()
                       if n.GetIdx() in ring_set]
        c1 = gly_n = None
        for c in ring_o_nbrs:
            for nbr in mol.GetAtomWithIdx(c).GetNeighbors():
                if nbr.GetIdx() not in ring_set and nbr.GetSymbol() == "N":
                    c1, gly_n = c, nbr.GetIdx()
                    break
            if c1 is not None:
                break
        if c1 is None:
            continue
        c4 = next(i for i in ring_o_nbrs if i != c1)
        c2 = next((n.GetIdx() for n in mol.GetAtomWithIdx(c1).GetNeighbors()
                    if n.GetIdx() in ring_set and n.GetIdx() != ring_o), None)
        c3 = next((n.GetIdx() for n in mol.GetAtomWithIdx(c4).GetNeighbors()
                    if n.GetIdx() in ring_set and n.GetIdx() != ring_o), None)
        c5 = next((n.GetIdx() for n in mol.GetAtomWithIdx(c4).GetNeighbors()
                    if n.GetIdx() not in ring_set and n.GetSymbol() == "C"), None)

        base = identify_thermo_base(mol, gly_n, ring_set)
        rib2 = identify_thermo_2prime(mol, c2, ring_set)
        base_atoms = _find_base_atoms(mol, gly_n, ring_set)

        p5 = _find_phosphorus_via_oxygen(mol, c5, ring_set) if c5 else None
        p3 = _find_phosphorus_via_oxygen(mol, c3, ring_set) if c3 else None

        nuc_atoms = ring_set | base_atoms
        if c5:
            nuc_atoms.add(c5)
        all_nuc_atoms |= nuc_atoms

        conj_2prime = detect_2prime_conjugate(mol, c2, ring_set)

        nucleotides.append({
            "base": base, "rib2": rib2,
            "p5": p5, "p3": p3,
            "c4": c4, "c5": c5,
            "ring_set": ring_set,
            "conj_2prime": conj_2prime,
        })

    # Order 5'->3'
    p5_map = {n["p5"]: i for i, n in enumerate(nucleotides) if n["p5"] is not None}
    all_p3 = {n["p3"] for n in nucleotides if n["p3"] is not None}
    all_p5 = {n["p5"] for n in nucleotides if n["p5"] is not None}
    start = next((i for i, n in enumerate(nucleotides) if n["p5"] is None), None)
    if start is None:
        start = next((i for i, n in enumerate(nucleotides) if n["p5"] not in all_p3), 0)
    ordered = []
    cur, seen = start, set()
    while cur is not None and cur not in seen:
        seen.add(cur)
        ordered.append(cur)
        cur = p5_map.get(nucleotides[cur]["p3"])

    # VP detection
    has_vp = False
    if ordered:
        first = nucleotides[ordered[0]]
        if first["p5"] is None:
            has_vp = detect_5prime_vp(mol, first)

    # Collect conjugates
    conjugates = []
    for pos, idx in enumerate(ordered):
        n = nucleotides[idx]
        if n["conj_2prime"] is not None:
            atoms, desc = n["conj_2prime"]
            formula = compute_fragment_formula(mol, atoms)
            cmono, cavg = compute_fragment_mw(mol, atoms)
            conjugates.append({
                "position": pos + 1, "location": "2'",
                "formula": formula, "description": desc,
                "mono_mw": cmono, "avg_mw": cavg
            })
    if ordered:
        last = nucleotides[ordered[-1]]
        if last["p3"] is not None:
            result = detect_3prime_conjugate(mol, last["p3"], all_nuc_atoms, all_p5)
            if result:
                atoms, desc = result
                formula = compute_fragment_formula(mol, atoms)
                cmono, cavg = compute_fragment_mw(mol, atoms)
                conjugates.append({
                    "position": len(ordered), "location": "3' terminal",
                    "formula": formula, "description": desc,
                    "mono_mw": cmono, "avg_mw": cavg
                })
        first = nucleotides[ordered[0]]
        result = detect_5prime_conjugate(mol, first, all_nuc_atoms)
        if result:
            atoms, desc = result
            formula = compute_fragment_formula(mol, atoms)
            cmono, cavg = compute_fragment_mw(mol, atoms)
            conjugates.append({
                "position": 1, "location": "5' terminal",
                "formula": formula, "description": desc,
                "mono_mw": cmono, "avg_mw": cavg
            })

    return {
        "mol": mol, "nucleotides": nucleotides, "ordered": ordered,
        "has_vp": has_vp, "sugar_like": sugar_like,
        "mono_mw": mono_mw, "avg_mw": avg_mw,
        "conjugates": conjugates, "n_nucleotides": len(ordered)
    }

# ============================================================
# Thermo FULL formatter
# ============================================================
def format_thermo(parsed):
    mol = parsed["mol"]
    nts = parsed["nucleotides"]
    ordered = parsed["ordered"]
    sugar_like = parsed["sugar_like"]
    if not ordered:
        return ""

    first = nts[ordered[0]]
    prefix = ""
    if parsed["has_vp"]:
        prefix = "v"
    elif first["p5"] is not None:
        prefix = classify_thermo_linker(mol, first["p5"], exclude=sugar_like) or ""

    out = []
    for pos, idx in enumerate(ordered):
        n = nts[idx]
        monomer = f"{n['base']}{n['rib2']}"
        if pos == 0:
            out.append(prefix + monomer)
        else:
            prev_p3 = nts[ordered[pos - 1]]["p3"]
            lk = classify_thermo_linker(mol, prev_p3, exclude=sugar_like) or "p"
            out.append(f"-{lk}{monomer}")

    thermo = "".join(out)

    # 3' terminal phosphate
    all_p5 = {n["p5"] for n in nts if n["p5"] is not None}
    last = nts[ordered[-1]]
    if last["p3"] is not None and last["p3"] not in all_p5:
        lk3 = classify_thermo_linker(mol, last["p3"], exclude=sugar_like)
        if lk3:
            thermo += f"-{lk3}"

    return thermo


# ============================================================
# Agilent formatter (UNCHANGED)
# ============================================================
def format_agilent(parsed):
    mol = parsed["mol"]
    nts = parsed["nucleotides"]
    ordered = parsed["ordered"]
    if not ordered:
        return ""

    parts = []
    first = nts[ordered[0]]
    if parsed["has_vp"]:
        parts.append("/5VP/")
    elif first["p5"] is not None:
        parts.append("/5Phos/")

    def agi_monomer(base, rib2):
        if rib2 == "m":
            return f"m{base}"
        if rib2 == "f":
            return f"f{base}"
        if rib2 == "r":
            return f"r{base}"
        if rib2 == "d":
            return base
        return base

    for pos, idx in enumerate(ordered):
        n = nts[idx]
        if pos > 0:
            prev_p3 = nts[ordered[pos - 1]]["p3"]
            if is_phosphorothioate(mol, prev_p3):
                parts.append("*")
        parts.append(agi_monomer(n["base"], n["rib2"]))

    return "".join(parts)


# ============================================================
# Top-level: parse + format everything
# ============================================================
def smiles_to_all(smiles_string):
    strand_smiles = smiles_string.split(".")
    n_strands = len(strand_smiles)

    # Full product MW
    full_mol = Chem.MolFromSmiles(smiles_string)
    if full_mol is None:
        raise ValueError("RDKit could not parse the SMILES string.")
    full_mono = Descriptors.ExactMolWt(full_mol)
    full_avg = Descriptors.MolWt(full_mol)

    # Parse each strand
    strand_data = []
    for s in strand_smiles:
        strand_data.append(parse_strand(s))

    # Sequences
    thermo_parts = [format_thermo(sd) for sd in strand_data]
    agilent_parts = [format_agilent(sd) for sd in strand_data]
    thermo_str = ".".join(thermo_parts)
    agilent_str = ".".join(agilent_parts)

    # MW summary
    mw_lines = []
    mw_lines.append("  Full-Length Product:")
    mw_lines.append(f"    Monoisotopic: {full_mono:.3f} Da")
    mw_lines.append(f"    Average:      {full_avg:.2f} Da")
    mw_lines.append("")
    for i, sd in enumerate(strand_data):
        label = f"Strand {i+1}"
        if n_strands == 2:
            label += " (sense)" if i == 0 else " (antisense)"
        mw_lines.append(f"  {label}  ({sd['n_nucleotides']} nt):")
        mw_lines.append(f"    Monoisotopic: {sd['mono_mw']:.3f} Da")
        mw_lines.append(f"    Average:      {sd['avg_mw']:.2f} Da")

    # Conjugates
    conj_lines = []
    all_conj = []
    for i, sd in enumerate(strand_data):
        for c in sd["conjugates"]:
            c["strand"] = i + 1
            all_conj.append(c)
    if all_conj:
        for c in all_conj:
            prefix = f"Strand {c['strand']}, " if n_strands > 1 else ""
            conj_lines.append(f"  {prefix}Position {c['position']} ({c['location']}): {c['formula']}")
            conj_lines.append(f"    Description:  {c['description']}")
            conj_lines.append(f"    Monoisotopic: {c['mono_mw']:.3f} Da")
            conj_lines.append(f"    Average:      {c['avg_mw']:.2f} Da")
    else:
        conj_lines.append("  None detected")

    return thermo_str, agilent_str, mw_lines, conj_lines


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    if len(sys.argv) > 1:
        SMILES = sys.argv[1]
    else:
        print("=" * 60)
        print("SMILES to Oligonucleotide Shorthand Converter")
        print("=" * 60)
        print("Paste SMILES (single or duplex with '.') then press Enter:")
        SMILES = input("SMILES> ").strip()

    if not SMILES:
        print("Error: No SMILES provided.")
        sys.exit(1)

    thermo, agilent, mw_lines, conj_lines = smiles_to_all(SMILES)

    print()
    print(f"Thermo BioPharma Finder Sequence: {thermo}")
    print(f"Agilent BioConfirm Sequence:      {agilent}")
    print()
    print("=" * 60)
    print("MOLECULAR WEIGHT SUMMARY")
    print("=" * 60)
    for line in mw_lines:
        print(line)
    print()
    print("=" * 60)
    print("CONJUGATE(S)")
    print("=" * 60)
    for line in conj_lines:
        print(line)