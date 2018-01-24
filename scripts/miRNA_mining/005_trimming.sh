#!/bin/bash
command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed. Aborting." >&2; exit 1; }

./cutadapt_folder.pl raw_smRNA/ trim_smRNA/