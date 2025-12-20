import subprocess
import argparse
from pathlib import Path
import sys
import pandas as pd

def run_flames_direct(args):
    """
    Run FLAMES using precomputed:
      - SuSiE credible set files (directory)
      - MAGMA gene-level output
      - MAGMA covariate output
      - PoPS score file
    """
    # ==============================================================
    # Helper: Create FLAMES indexfile
    # ==============================================================
    def create_flames_indexfile(finemap_dir: Path, outdir: Path, indexfile_path: str):
        finemap_dir = Path(finemap_dir)

        if not finemap_dir.exists() or not finemap_dir.is_dir():
            raise NotADirectoryError(
                f"ERROR: --finemap_cred_dir must be a folder, not a file.\n"
                f"Provided: {finemap_dir}"
            )

        outdir.mkdir(parents=True, exist_ok=True)

        # Collect SuSiE credible set files
        finemap_files = sorted([
            f for f in finemap_dir.iterdir()
            if f.is_file() and not f.name.startswith(".")
        ])

        if len(finemap_files) == 0:
            raise RuntimeError(
                f"No SuSiE credible-set files found in folder: {finemap_dir}"
            )

        rows = []
        for f in finemap_files:
            prefix = f.stem
            annotfile = outdir / f"{prefix}_locus_1.txt"
            rows.append({
                "Filename": str(f.resolve()),
                "Annotfiles": str(annotfile.resolve())
            })

        df = pd.DataFrame(rows)
        df.to_csv(indexfile_path, sep="\t", index=False)

        return indexfile_path

    # ==============================================================
    # PREPARE DIRECT INPUTS
    # ==============================================================

    # Folder containing ALL SuSiE credible-set files
    finemap_dir = Path(args.finemap_cred_dir)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Create indexfile
    indexfile_path = outdir / f"{args.sample_id}_flames_index.tsv"

    indexfile = create_flames_indexfile(
        finemap_dir=finemap_dir,
        outdir=outdir,
        indexfile_path=str(indexfile_path)
    )

    # FLAMES inputs
    annotation_dir = args.flames_annot_dir
    pops          = args.pops_score_file
    magma_z       = args.magma_genes_out
    magma_tissue  = args.magma_tissue_covar_results

    # Defaults for FLAMES model
    filename  = "FLAMES_scores"
    modelpath = Path(__file__).parent / "model"
    snp_col   = "cred1"
    prob_col  = "prob1"
    distance  = 750000
    weight    = 0.725

    flames_py = Path(__file__).parent / "FLAMES.py"

    # ==============================================================
    # STEP 1 — FLAMES annotate
    # ==============================================================
    cmd_annot = f"""
    python "{flames_py}" annotate \
        --annotation_dir "{annotation_dir}" \
        --pops "{pops}" \
        --magma_z "{magma_z}" \
        --magma_tissue "{magma_tissue}" \
        --indexfile "{indexfile}" \
        --prob_col "{prob_col}" \
        --SNP_col "{snp_col}"
    """

    print("             🔥 Running FLAMES annotation…")
    subprocess.run(cmd_annot, shell=True, check=True)

    # ==============================================================
    # STEP 2 — FLAMES scoring
    # ==============================================================
    cmd_score = f"""
    python "{flames_py}" FLAMES \
        --indexfile "{indexfile}" \
        --outdir "{outdir}" \
        --filename "{filename}" \
        --distance {distance} \
        --weight {weight} \
        --modelpath "{modelpath}"
    """
    print("             🔥 Running FLAMES scoring…")
    subprocess.run(cmd_score, shell=True, check=True)
    print("             🎉 FLAMES Analysis completed successfully!")
    print(" ")

