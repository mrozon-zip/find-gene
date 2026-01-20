#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List

OUTFMT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen"
]

def run(cmd: List[str], cwd: str | None = None) -> None:
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

def shell_out(cmd: List[str]) -> str:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout

def split_commas(values: List[str]) -> List[str]:
    out: List[str] = []
    for v in values:
        out.extend([p.strip() for p in v.split(",") if p.strip()])
    return out

def read_list_files(files: List[str]) -> List[str]:
    items: List[str] = []
    for fp in files:
        p = Path(fp)
        if not p.exists():
            raise FileNotFoundError(f"List file not found: {fp}")
        for line in p.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            items.append(s)
    return items

def sanitize_basename(path_str: str) -> str:
    name = Path(path_str).name
    name = re.sub(r"\.(fa|fna|fasta|fas)(\.gz)?$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name

def ensure_fasta_download(get_ref_path: str, accession: str, queries_dir: Path) -> Path:
    queries_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = queries_dir / f"{accession}.fasta"

    # Cache hit if exists and non-empty
    if fasta_path.exists() and fasta_path.stat().st_size > 0:
        return fasta_path

    fasta = shell_out([get_ref_path, accession])
    if not fasta.startswith(">"):
        raise RuntimeError(
            f"Downloaded content for {accession} does not look like FASTA.\n"
            f"First 200 chars:\n{fasta[:200]}"
        )
    fasta_path.write_text(fasta)
    return fasta_path

def db_exists(db_prefix: Path) -> bool:
    # For nucleotide BLAST DB, core files usually include:
    # .nhr .nin .nsq (or sometimes .ndb etc depending on version)
    needed = [db_prefix.with_suffix(suf) for suf in [".nhr", ".nin", ".nsq"]]
    return all(p.exists() and p.stat().st_size > 0 for p in needed)

def make_db_if_needed(assembly_path: Path, db_prefix: Path) -> None:
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    if db_exists(db_prefix):
        return  # cache hit

    run([
        "makeblastdb",
        "-in", str(assembly_path),
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", str(db_prefix)
    ])

def blast_one(query_fasta: Path, db_prefix: Path, out_tsv: Path) -> None:
    outfmt = "6 " + " ".join(OUTFMT_COLUMNS)
    run([
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-outfmt", outfmt,
        "-max_target_seqs", "20",
        "-evalue", "1e-20",
        "-out", str(out_tsv)
    ])

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Download NCBI accession FASTA, cache BLAST DB(s) from assembly fasta(s), run blastn, write cumulative TSV."
    )

    ap.add_argument("-a", "--assembly", action="append", default=[],
                    help="Assembly FASTA path. Repeatable. Comma-separated also accepted.")
    ap.add_argument("-A", "--assembly-file", action="append", default=[],
                    help="Text file listing assembly FASTA paths, one per line (blank/# ignored). Repeatable.")

    ap.add_argument("-i", "--accession", action="append", default=[],
                    help="NCBI accession. Repeatable. Comma-separated also accepted.")
    ap.add_argument("-I", "--accession-file", action="append", default=[],
                    help="Text file listing accessions, one per line (blank/# ignored). Repeatable.")

    ap.add_argument("--get-ref", default="/usr/local/bin/get_ref.sh",
                    help="Path to get_ref.sh inside container.")
    ap.add_argument("--workdir", default=".work",
                    help="Workdir for cached queries and BLAST DBs (default: ./.work).")
    ap.add_argument("--outdir", default="output",
                    help="Output directory in current working directory (default: ./output).")
    ap.add_argument("-o", "--output", default="cumulative_results.tsv",
                    help="Output TSV filename (written inside outdir).")

    args = ap.parse_args()

    assemblies = split_commas(args.assembly) + read_list_files(args.assembly_file)
    accessions = split_commas(args.accession) + read_list_files(args.accession_file)

    if not assemblies:
        raise ValueError("No assemblies provided. Use -a and/or -A.")
    if not accessions:
        raise ValueError("No accessions provided. Use -i and/or -I.")

    # Resolve paths
    assemblies_paths = [Path(a).resolve() for a in assemblies]
    for p in assemblies_paths:
        if not p.exists():
            raise FileNotFoundError(f"Assembly not found: {p}")

    workdir = Path(args.workdir).resolve()
    queries_dir = workdir / "queries"
    db_dir = workdir / "blastdb"
    tmp_dir = workdir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / args.output

    # Download/cached queries
    accession_to_fasta: dict[str, Path] = {}
    for acc in accessions:
        accession_to_fasta[acc] = ensure_fasta_download(args.get_ref, acc, queries_dir)

    # Build/cached DBs
    assembly_to_dbprefix: dict[Path, Path] = {}
    for asm_path in assemblies_paths:
        asm_name = sanitize_basename(str(asm_path))
        db_prefix = db_dir / f"{asm_name}_asm_db"
        make_db_if_needed(asm_path, db_prefix)
        assembly_to_dbprefix[asm_path] = db_prefix

    # Write cumulative TSV
    header = ["comparison"] + OUTFMT_COLUMNS

    with out_path.open("w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(header)

        for acc in accessions:
            qf = accession_to_fasta[acc]
            for asm_path in assemblies_paths:
                asm_name = sanitize_basename(str(asm_path))
                comparison = f"{acc}_vs_{asm_name}_asm_db"

                db_prefix = assembly_to_dbprefix[asm_path]
                tmp_out = tmp_dir / f"{comparison}.tsv"

                blast_one(qf, db_prefix, tmp_out)

                if tmp_out.exists() and tmp_out.stat().st_size > 0:
                    for line in tmp_out.read_text().splitlines():
                        if not line.strip():
                            continue
                        writer.writerow([comparison] + line.split("\t"))
                else:
                    # Keep one row per comparison even if no hits
                    writer.writerow([comparison] + [""] * len(OUTFMT_COLUMNS))

    print(f"Done. Cumulative results: {out_path}")
    print(f"Cached workdir: {workdir}")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
