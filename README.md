Here is a description of a workflow that can be used for determination if a specific sequence is found in a database (database in this context is an assembly of a genome)

- [Requirements:](#requirements)
- [Environment:](#environment)
- [Workflow steps](#workflow-steps)
- [Results description](#results-description)
- [How to draw a biological conclusion](#how-to-draw-a-biological-conclusion)
	- [Example interpretation sentence](#example-interpretation-sentence)
- [Usage](#usage)
- [Inputs](#inputs)
- [Multiple inputs](#multiple-inputs)
- [Output](#output)
	- [Output directory](#output-directory)
	- [Output file format (`cumulative_results.tsv`)](#output-file-format-cumulative_resultstsv)
- [Caching behavior (important)](#caching-behavior-important)
- [Help](#help)

# Workflow
## Requirements:
1. Assembled genome (fasta format)
2. Accession number for a query sequence
3. BLAST environment
***
## Environment:
```
name: blast_env
channels:
  - defaults
  - bioconda
  - conda-forge
dependencies:
  - _libgcc_mutex=0.1
  - _openmp_mutex=5.1
  - blast=2.17.0
  - bzip2=1.0.8
  - c-ares=1.34.5
  - ca-certificates=2025.12.2
  - curl=8.17.0
  - entrez-direct=24.0
  - gettext=0.21.0
  - icu=73.1
  - jansson=2.14
  - libcurl=8.17.0
  - libev=4.33
  - libgcc=15.2.0
  - libgcc-ng=15.2.0
  - libgomp=15.2.0
  - libiconv=1.16
  - libidn2=2.3.8
  - libnghttp2=1.67.1
  - libnsl=2.0.0
  - libsqlite=3.51.2
  - libssh2=1.11.1
  - libstdcxx=15.2.0
  - libstdcxx-ng=15.2.0
  - libunistring=1.3
  - libxml2=2.13.9
  - libzlib=1.3.1
  - lz4-c=1.9.4
  - ncbi-vdb=3.3.0
  - ncurses=6.5
  - openssl=3.0.18
  - pcre2=10.46
  - perl=5.32.1
  - perl-archive-tar=3.04
  - perl-carp=1.38
  - perl-common-sense=3.75
  - perl-compress-raw-bzip2=2.214
  - perl-compress-raw-zlib=2.214
  - perl-encode=3.19
  - perl-exporter=5.72
  - perl-exporter-tiny=1.002002
  - perl-extutils-makemaker=7.70
  - perl-io-compress=2.214
  - perl-io-zlib=1.15
  - perl-json=4.10
  - perl-json-xs=4.04
  - perl-list-moreutils=0.430
  - perl-list-moreutils-xs=0.430
  - perl-parent=0.236
  - perl-pathtools=3.75
  - perl-scalar-list-utils=1.62
  - perl-types-serialiser=1.01
  - wget=1.25.0
  - xz=5.6.4
  - zlib=1.3.1
  - zstd=1.5.7
```

Constructed using:
```sh
conda create -n blast_env -c bioconda blast
conda activate blast_env 
```

## Workflow steps
1. Create DB

```
makeblastdb -in <assembly_input.fasta> -dbtype nucl -parse_seqids -out <db_name>
```

2. Query

```
blastn -query <query.fasta> -db <db_name>   -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen qcovs"   -max_target_seqs 5 -evalue 1e-20   > <results.tsv>
```
***
## Results description
The file `results.tsv` contains the results of a BLASTN search, where a **reference nucleotide sequence** (query) was compared against an **assembled bacterial genome** (subject).  
Each row corresponds to **one alignment (hit)** between the query sequence and a region of the assembled genome.

1. `qseqid` — Query sequence ID
This indicates _which reference sequence_ was searched for in the genome.

2. `sseqid` — Subject sequence ID
This tells you _where in the assembly_ the matching sequence is located.  
Different `sseqid` values correspond to different contigs (chromosome fragments, plasmids, or unplaced contigs).

3. `pident` — Percent identity
- **~100%** → nearly identical sequence    
- **95–99%** → very close homolog (same gene, different strain)
- **<90%** → more distant homolog    
This value reflects **sequence similarity**, not gene expression or functionality.

 4. `length` — Alignment length (bp)
Long alignments covering most or all of the query sequence provide strong evidence that the gene or fragment is present.

5. `qstart` — Query start position
Position in the query sequence where the alignment begins.

6. `qend` — Query end position
These columns show **which part of the reference sequence** aligns to the genome.  
If `qstart = 1` and `qend ≈ query length`, the **entire reference sequence is covered**.

7. `sstart` — Subject start position
Genomic coordinate (on the contig) where the alignment begins.

8. `send` — Subject end position
These columns define the **exact genomic location** of the matching sequence in the assembly and can be used to extract the sequence or inspect its genomic context.

9. `evalue` — Expect value
- **0.0 or extremely small values (e.g. 1e-100)** → match is not due to chance
- Larger values → less reliable matches
For bacterial gene searches, an e-value of **0.0** indicates a **highly significant hit**

10. `bitscore` — Bit score
Higher bit scores indicate better alignments.  
Bit scores are mainly useful for **ranking hits**, not biological interpretation on their own.

11. `qlen` — Query length (bp) _(if present)_
Used to assess how much of the reference sequence is represented in the alignment.

12. `qcovs` — Query coverage (%) _(if present)_
- **~100%** → full-length reference sequence is present
- **<100%** → only a fragment is present
This is crucial for deciding whether a **full gene** or only a **partial sequence** is detected.
***
## How to draw a biological conclusion

A **strong conclusion that the reference sequence is present in the genome** can be made when:
- Percent identity is high (typically >95%)
- Query coverage is high (ideally ~100%)
- Alignment length matches the expected gene length
- E-value is extremely small (≈ 0)

A **weaker but still meaningful conclusion** (presence of a homolog) applies when:
- Identity is moderately high
- Coverage is partial
- The organism is evolutionarily related

Importantly:
	This analysis demonstrates **sequence presence**, not gene expression, regulation, or functionality.
***
### Example interpretation sentence 
"BLASTN analysis revealed a high-confidence match between the reference sequence and contig_1 of the assembled genome, with 98.8% nucleotide identity over a 2,535 bp alignment and full query coverage, indicating the presence of a closely related homolog of the reference gene in the analyzed bacterial genome."

# Container

`find-gene.sif` is an Apptainer container that checks whether one or more **reference nucleotide sequences** (given by NCBI accession IDs) are present in one or more **assembled bacterial genomes**.

The container:

- downloads reference sequences from NCBI (cached locally),
- builds BLAST nucleotide databases from provided assemblies (cached locally),
- runs `blastn` for all accession × assembly combinations,
- outputs **one cumulative TSV file** with all results

The container answers the question:
> _“Is sequence A present (or highly similar) in genome B?”_

## Usage

Simple usage:
```
apptainer run find-gene.sif \
  -a genome.fasta \
  -i EU188855.1
```

Output:
```
./output/cumulative_results.tsv
```

## Inputs
Both inputs are required

Assemblies - must be fasta file format representing assembled genomes or contigs
Three posibilities on how to provide them:
1. Repeating `-a` flag (recommended for small numbers)
```
-a asm1.fasta -a asm2.fasta
```
2. Assembly list file (recommended for many assemblies)
Create `assemblies.txt`:
```
asm1.fasta
asm2.fasta
```
Use:
```
-A assemblies.txt
```
3. Comma-separated (supported, less readable)
```
-a asm1.fasta,asm2.fasta
```

Accession IDs - Accession IDs must be valid **NCBI nucleotide accessions**.

can be provided in a same way as assemblies:
```
-i EU188855.1 -i JX398137.1
```

```
-I accessions.txt
```

with accessions.txt:
```
JX398137.1
EU188855.1
```

```
-i EU188855.1,JX398137.1
```

## Multiple inputs
The container automatically performs a **Cartesian product**:
> _N accessions × M assemblies = N×M BLAST comparisons_

## Output

### Output directory

All results are written to `./output/` Default output file:`output/cumulative_results.tsv` You can rename it: `-o my_results.tsv`

### Output file format (`cumulative_results.tsv`)

Tab-separated table with a header row.

| Column       | Meaning                                                       |
| ------------ | ------------------------------------------------------------- |
| `comparison` | Identifier of the comparison (`ACCESSION_vs_ASSEMBLY_asm_db`) |
| `qseqid`     | Query accession ID                                            |
| `sseqid`     | Contig ID in the assembly                                     |
| `pident`     | Percent nucleotide identity                                   |
| `length`     | Alignment length (bp)                                         |
| `qstart`     | Alignment start in query                                      |
| `qend`       | Alignment end in query                                        |
| `sstart`     | Alignment start in assembly                                   |
| `send`       | Alignment end in assembly                                     |
| `evalue`     | Statistical significance                                      |
| `bitscore`   | Alignment score                                               |
| `qlen`       | Query sequence length                                         |

Each row represents **one BLAST hit**. If no hit is found for a given comparison, a row with empty alignment fields is still written (so all combinations are represented).

## Caching behavior (important)

The container creates a hidden working directory: `./.work/`

Inside:

```.work/
├── queries/     # downloaded accession FASTA files   
└── blastdb/     # BLAST databases built from assemblies```
```

This means:
- Reference sequences are downloaded **once**
- BLAST databases are built **once per assembly**
- Subsequent runs are much faster
You may delete `.work/` at any time to force a clean re-run.
## Help
To see all options:
`apptainer run find-gene.sif --help`
