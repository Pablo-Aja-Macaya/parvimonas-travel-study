import pandas as pd
import os
from subprocess import run

os.chdir("path/to/Repository/genomics/comparative_genomics/crisprs")

# Merge crosses in data/cross into one file
cmd = "cat data/cross/PM89KC_AC_1_39528__cross_*k__prot.fasta > data/cross/merged.fasta"
run(cmd, shell=True, check=True)

# Make DIAMOND database with merged file
db_fasta = "data/cross/merged.fasta"
database = "data/cross/database.dmnd"
outfmt_fields = "qseqid qtitle sseqid stitle pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore"
diamond_program = "blastp"
outfmt = 6
max_target_seqs = 10000000000000
threads = 60
max_hsps = 1

cmd = f"diamond makedb  --in {db_fasta} --db {database}"
print(cmd)
run(cmd, shell=True, check=True)

# Query the merged file against itself
query = "data/cross/merged.fasta"
output_file = "data/cross/dmnd.tsv"
cmd = f"""
diamond {diamond_program} --in {db_fasta} --db {database} \
    --query {query} --outfmt {outfmt} {outfmt_fields} \
    -o {output_file} --max-target-seqs {max_target_seqs} \
    --max-hsps {max_hsps} --threads {threads} --salltitles \
    --quiet --more-sensitive
"""
print(cmd)
run(cmd, shell=True, check=True)

# Read results
dmnd_df = pd.read_csv(output_file, sep="\t", header=None)
dmnd_df.columns = outfmt_fields.split()

# Remove self hits
dmnd_df = dmnd_df[dmnd_df["qseqid"] != dmnd_df["sseqid"]]
dmnd_df.sort_values(["qseqid"])

# Save results
dmnd_df.to_csv("data/cross/dmnd_proc.tsv", sep="\t", index=None)
