#!/usr/bin/env bash

set -euo pipefail

# ---- defaults ----
DB="nuccore"
RETTYPE="fasta"
RETMODE="text"
BASE_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ---- argument check ----
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accession_id>" >&2
    exit 1
fi

ID="$1"

# ---- fetch ----
curl -L \
  "${BASE_URL}?db=${DB}&id=${ID}&rettype=${RETTYPE}&retmode=${RETMODE}"
