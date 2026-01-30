"""
Analyse all PDB files in the pdb folder: number of chains, total amino acids,
experiment type (NMR vs single/model), and write results to CSV.
"""
import os
import csv
from pathlib import Path

# Standard 20 amino acid residue names (3-letter)
STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})


def analyse_pdb(filepath):
    """
    Parse a PDB file and return num_chains, total_amino_acids, expdta, num_models.
    Only counts standard amino acids from ATOM records (first model if multi-model).
    """
    residues_seen = set()  # (chain_id, res_seq) for standard AA
    expdta = ""
    num_models = 0
    model_count = 0  # number of MODEL lines seen
    count_atoms = True  # count residues only in first model

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("EXPDTA"):
                expdta = line[10:80].strip()
            elif line.startswith("NUMMDL"):
                try:
                    num_models = int(line[10:14].strip())
                except ValueError:
                    pass
            elif line.startswith("MODEL"):
                model_count += 1
                if num_models == 0:
                    num_models = 1
                # Stop counting residues after first model (when we see 2nd MODEL)
                if model_count > 1:
                    count_atoms = False
            elif line.startswith("ENDMDL"):
                count_atoms = False  # only count first model's ATOMs
            elif line.startswith("ATOM") and count_atoms:
                res_name = line[17:20].strip()
                chain_id = line[21]
                try:
                    res_seq = int(line[22:26].strip())
                except ValueError:
                    continue
                if res_name in STANDARD_AA:
                    residues_seen.add((chain_id, res_seq))

    # Use MODEL count if we saw MODEL lines; else use NUMMDL or 1
    if model_count > 0:
        num_models = model_count
    elif num_models == 0:
        num_models = 1  # single structure, no MODEL records

    chains = {c for c, _ in residues_seen}
    chain_names = ",".join(sorted(chains)) if chains else ""
    num_chains = len(chains)
    total_amino_acids = len(residues_seen)

    # Experiment type: NMR vs single (X-ray, EM, etc.)
    exp_upper = expdta.upper()
    if "NMR" in exp_upper or "NUCLEAR MAGNETIC" in exp_upper:
        structure_type = "NMR"
    elif num_models > 1:
        structure_type = "NMR"  # multi-model ensemble typically NMR
    else:
        structure_type = "single"

    return {
        "pdb_id": Path(filepath).stem.lower(),
        "num_chains": num_chains,
        "chain_names": chain_names,
        "total_amino_acids": total_amino_acids,
        "structure_type": structure_type,
        "expdta": expdta,
        "num_models": num_models,
    }


def run_analysis(pdb_dir="pdb", output_csv="pdb_analysis.csv"):
    """Analyse all .pdb files in pdb_dir and write CSV."""
    pdb_dir = Path(pdb_dir)
    if not pdb_dir.is_dir():
        raise FileNotFoundError(f"Not a directory: {pdb_dir}")

    results = []
    for path in sorted(pdb_dir.glob("*.pdb")):
        try:
            row = analyse_pdb(path)
            results.append(row)
        except Exception as e:
            results.append({
                "pdb_id": path.stem.lower(),
                "num_chains": "",
                "chain_names": "",
                "total_amino_acids": "",
                "structure_type": "error",
                "expdta": str(e),
                "num_models": "",
            })

    fieldnames = ["pdb_id", "num_chains", "chain_names", "total_amino_acids", "structure_type", "expdta", "num_models"]
    outpath = Path(output_csv)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)

    return results, str(outpath)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyse PDBs: chains, amino acids, NMR vs single")
    parser.add_argument("--pdb-dir", default="pdb", help="Directory containing .pdb files")
    parser.add_argument("--output", "-o", default="pdb_analysis.csv", help="Output CSV path")
    args = parser.parse_args()
    results, outpath = run_analysis(pdb_dir=args.pdb_dir, output_csv=args.output)
    print(f"Analysed {len(results)} PDB files. Results written to {outpath}")
