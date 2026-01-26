
# -*- coding: utf-8 -*-
"""
Run LDSC (HM3 filter → munge → h2) on Windows WITHOUT --merge-alleles

Pipeline:
  Step 0) Filter GWAS summary statistics to HapMap3 rsIDs (chunked)
  Step 1) munge_sumstats.py
          - --no-alleles
          - --ignore BETA,SE,Z
          - --a1-inc
  Step 2) ldsc.py --h2

Rationale:
  This pipeline avoids allele merging and strand alignment at the LDSC stage.
  SNP harmonization (rsID, REF/ALT, effect direction) is assumed to be handled
  upstream during GWAS QC and summary-stat preparation.

Author : Tung-Yen Wu
"""

import subprocess
import sys
import gzip
from pathlib import Path
import pandas as pd


# ============================================================
# 0) USER CONFIG — EDIT THIS SECTION ONLY
# ============================================================

CFG = dict(

    # ---- LDSC installation ----
    LDSC_DIR = Path(r"C:\PATH\TO\ldsc"),
    LDSC_PY  = "ldsc.py",
    MUNGE_PY = "munge_sumstats.py",

    # ---- Input GWAS summary statistics ----
    SUMSTATS = Path(r"C:\PATH\TO\resilience_ldsc_ready.sumstats.gz"),

    # ---- HapMap3 SNP list ----
    # Can be:
    #   - rsID-only list (hapmap_3.3.hg38.txt)
    #   - SNP A1 A2 format (hapmap_3.3.hg38.alleles.txt)
    HM3_LIST = Path(r"C:\PATH\TO\hm3.snplist"),

    # ---- LD score reference (hg38, 1000G Phase 3) ----
    LD38_DIR = Path(
        r"C:\PATH\TO\LDSCORE_1000G_Phase3_ldscores\LDscore\LDscore."
    ),

    # ---- Output directory ----
    OUTDIR = Path(r"C:\PATH\TO\ldsc_out"),

    # ---- Column names in sumstats ----
    SNP_COL = "SNP",
    A1_COL  = "A1",
    A2_COL  = "A2",
    P_COL   = "P",
    N_COL   = "N",

    # ---- Runtime parameters ----
    CHUNK_SIZE = 500_000
)


# ============================================================
# 1) Utility functions
# ============================================================

PYTHON = sys.executable


def run(cmd, cwd=None):
    """Run external command and stream output."""
    print(">>", " ".join(str(x) for x in cmd))
    proc = subprocess.Popen(
        cmd, cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    for line in proc.stdout:
        print(line, end="")
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"Command failed (exit {ret})")


def preview_header(path, n=1):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            print("[Preview]", line.strip())


def load_hm3_rsids(hm3_path):
    """
    Load HapMap3 rsID set.
    Supports:
      - SNP
      - ID
      - VCF-like formats (3rd column as ID)
    """
    df = pd.read_csv(hm3_path, sep=r"\s+")
    cols = [c.upper() for c in df.columns]

    if "SNP" in cols:
        rs = df.iloc[:, cols.index("SNP")]
    elif "ID" in cols:
        rs = df.iloc[:, cols.index("ID")]
    else:
        rs = df.iloc[:, 2]  # fallback for VCF-like format

    rs = rs.astype(str)
    rs = rs[rs.str.startswith("rs")]
    return set(rs.values)


def filter_sumstats_to_hm3(sumstats, hm3_rs, out_path, snp_col, chunk_size):
    """
    Chunk-wise filtering of sumstats to HM3 rsIDs.
    """
    opener = gzip.open if sumstats.suffix == ".gz" else open

    with opener(sumstats, "rt", encoding="utf-8") as f:
        header = f.readline().strip().split()

    if snp_col not in header:
        raise ValueError(f"Column {snp_col} not found in sumstats")

    it = pd.read_csv(
        sumstats,
        sep=r"\s+",
        compression="infer",
        chunksize=chunk_size
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    first = True

    with gzip.open(out_path, "wt", encoding="utf-8") as g:
        for chunk in it:
            keep = chunk[chunk[snp_col].astype(str).isin(hm3_rs)]
            if keep.empty:
                continue
            keep.to_csv(
                g,
                sep="\t",
                index=False,
                header=first
            )
            first = False

    if first:
        raise RuntimeError("HM3 intersection is empty. Check rsID format and genome build.")

    return out_path


# ============================================================
# 2) Build LDSC commands
# ============================================================

def build_munge_cmd(in_sumstats, out_prefix):
    return [
        PYTHON, str(CFG["LDSC_DIR"] / CFG["MUNGE_PY"]),
        "--sumstats", str(in_sumstats),
        "--out", str(out_prefix),
        "--chunksize", str(CFG["CHUNK_SIZE"]),
        "--snp", CFG["SNP_COL"],
        "--a1",  CFG["A1_COL"],
        "--a2",  CFG["A2_COL"],
        "--p",   CFG["P_COL"],
        "--N-col", CFG["N_COL"],
        "--ignore", "BETA,SE,Z",
        "--a1-inc",
        "--no-alleles"
    ]


def build_h2_cmd(munged_sumstats, out_prefix):
    return [
        PYTHON, str(CFG["LDSC_DIR"] / CFG["LDSC_PY"]),
        "--h2", str(munged_sumstats),
        "--ref-ld-chr", str(CFG["LD38_DIR"]),
        "--w-ld-chr",   str(CFG["LD38_DIR"]),
        "--out", str(out_prefix)
    ]


# ============================================================
# 3) Main pipeline
# ============================================================

def main():

    print("\n== Preview summary statistics header ==")
    preview_header(CFG["SUMSTATS"])

    # Step 0: HM3 filter
    print("\n== Step 0: Filter to HapMap3 rsIDs ==")
    hm3_rs = load_hm3_rsids(CFG["HM3_LIST"])
    hm3_only = CFG["OUTDIR"] / "trait.hm3only.sumstats.gz"

    filter_sumstats_to_hm3(
        CFG["SUMSTATS"],
        hm3_rs,
        hm3_only,
        CFG["SNP_COL"],
        CFG["CHUNK_SIZE"]
    )

    print(f"[ok] HM3-only sumstats written: {hm3_only}")

    # Step 1: munge
    print("\n== Step 1: munge_sumstats.py ==")
    munged_prefix = CFG["OUTDIR"] / "trait.hm3.hg38"
    run(build_munge_cmd(hm3_only, munged_prefix))

    munged_sumstats = Path(str(munged_prefix) + ".sumstats.gz")
    if not munged_sumstats.exists():
        raise FileNotFoundError("munge output not found")

    # Step 2: LDSC h2
    print("\n== Step 2: LDSC h2 estimation ==")
    h2_prefix = CFG["OUTDIR"] / "trait_h2_hg38_1000G"
    run(build_h2_cmd(munged_sumstats, h2_prefix))

    print("\n== DONE ==")
    print(f"HM3-only : {hm3_only}")
    print(f"Munged   : {munged_sumstats}")
    print(f"LDSC out : {h2_prefix}.log / .results")


if __name__ == "__main__":
    main()