import pandas as pd
from subprocess import run
import os

os.chdir("path/to/Repository/genomics/comparative_genomics/crisprs")

############################################################
# ---------------------- Params -------------------------- #
############################################################
OUTFMT = "6 qseqid sseqid stitle sstrand pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qseq sseq"

# DB files
gpd_meta = "data/dbs/gpd/GPD_metadata.tsv"
gpd_seqs = "data/dbs/gpd/GPD_sequences.fa"
ovd_meta = "data/dbs/ovd/OVD-info.xlsx"
ovd_seqs = "data/dbs/ovd/OVD-genomes.fa"
merged_blast_db = "data/dbs/gpd_ovd/ovd__gpd.fasta"

# Queries
merged_spacers = "data/results/merged_spacers.fasta"
left_phage = "data/cross/PM89KC_AC_1__left_cross.fa"

# Output
out_folder = "output"

# Extra requirements: BLASTn in path


############################################################
# -------------------- Functions ------------------------- #
############################################################
def get_blast_results(blast_result, outfmt):
    df = pd.read_csv(
        blast_result,
        sep="\t",
        header=None,
    )
    df.columns = [i for i in outfmt.split(" ")[1:]]

    # Coverage
    df["cov"] = 100 * df["length"] / df["qlen"]

    # Sort
    df = df.sort_values(by=["sseqid", "cov", "pident"], ascending=False)
    df = df[["qseqid", "sseqid", "pident", "cov", "length", "qlen"]]

    return df


def get_last_host_taxon(row, sep, col):
    host_range_taxon = row[col]
    last = None
    if type(host_range_taxon) == str:
        host_range_taxon = host_range_taxon.split(sep)
        last = "; ".join(host_range_taxon[-2:])
    row["last_host_taxon"] = last
    return row


def merge_dbs(gpd_main_df, ovd_main_df, ovd_host_pred_df):
    """
    Merge databases GPD and OVD
    """
    # Process gpd metadata
    processed_gpd = gpd_main_df[["GPD_id", "Predicted_phage_taxon", "Host_range_taxon"]]
    processed_gpd.columns = ["sseqid", "phage_taxon", "host_range_taxon"]
    processed_gpd["db"] = "intestinal"

    # Process OVD metadata
    processed_ovd = pd.merge(ovd_main_df, ovd_host_pred_df, on="vOTU ID")
    processed_ovd = processed_ovd[
        [
            "vOTU ID",
            "Family-level taxonomy",
            "Host_taxonomy",
        ]
    ]
    processed_ovd.columns = ["sseqid", "phage_taxon", "host_range_taxon"]
    processed_ovd["db"] = "oral"

    # Merge metadata
    merged_ovd_gpd = pd.concat([processed_gpd, processed_ovd])
    merged_ovd_gpd["last_host_taxon"] = None
    merged_ovd_gpd["host_range_taxon"] = merged_ovd_gpd["host_range_taxon"].str.replace(
        "/", ";"
    )
    merged_ovd_gpd = merged_ovd_gpd.apply(
        get_last_host_taxon, sep=";", col="host_range_taxon", axis=1
    )
    merged_ovd_gpd = merged_ovd_gpd[["sseqid", "phage_taxon", "last_host_taxon", "db"]]

    return processed_gpd, processed_ovd, merged_ovd_gpd


############################################################
# ------------------ Set up metadata --------------------- #
############################################################

# Read OVD
ovd_meta_HP_df = pd.read_excel(ovd_meta, sheet_name="Host prediction")
ovd_meta_main_df = pd.read_excel(ovd_meta, sheet_name="Sheet1", skiprows=3, header=None)
ovd_meta_main_df.columns = [
    "vOTU ID",
    "Family-level taxonomy",
    "Family-level clade",
    "Genus-level clade",
    "Length (bp)",
    "Number of genes",
    "n_viral_genes__checkv",
    "n_microbial_genes__checkv",
    "quality__checkv",
    "completeness__checkv",
    "completeness_method__checkv",
    "contamination__checkv",
    "lytic/lysogenic__vibrant",
    "score__depvirfinder",
    "pvalue__depvirfinder",
    "max_score__virsorter2",
    "max_score_group__virsorter2",
]

# Read GPD
gpd_meta_df = pd.read_csv(gpd_meta, sep="\t")
gpd_meta_df = gpd_meta_df[
    ["GPD_id", "Source", "Size", "Predicted_phage_taxon", "Host_range_taxon"]
]

############################################################
# ----------- Make blast DB of OVD and GPD --------------- #
############################################################
cmd = f"makeblastdb -in {merged_blast_db} -dbtype nucl"
print(cmd)
run(cmd, shell=True, check=True)

############################################################
# --------- Blast spacers against OVD and GPD ------------ #
############################################################

# Run blast query
query = merged_spacers
db = merged_blast_db
out = f"{out_folder}/spacers_against_db.tsv"  # "/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/manual/KellyParvimonas_crispr/crispr/spacer_to_db/res.tsv"
outfmt = OUTFMT

cmd = [
    "blastn",
    f"-query {query}",
    f"-out {out}",
    f"-db {db}",
    "-perc_identity 90 -max_hsps 1 -max_target_seqs 1000 -num_threads 60 -task blastn-short",
    f"-outfmt '{outfmt}'",
]
cmd = " ".join(cmd)
print(cmd)
run(cmd, shell=True, check=True)

# Get blast results
merged_res = get_blast_results(out, OUTFMT)
merged_res = merged_res[merged_res["length"] >= merged_res["qlen"] * 0.8]

# Get merged database metadata
processed_gpd, processed_ovd, merged_ovd_gpd = merge_dbs(
    gpd_main_df=gpd_meta_df,
    ovd_main_df=ovd_meta_main_df,
    ovd_host_pred_df=ovd_meta_HP_df,
)

# Merge results with metadata
df = pd.merge(merged_res, merged_ovd_gpd, on="sseqid")
df[["qseqid_bc", "qseqid_crispr"]] = df["qseqid"].str.split("__", expand=True)
df = df.fillna("UNKNOWN")
df = (
    df.groupby(
        by=[
            "qseqid_bc",
            "qseqid_crispr",
            "sseqid",
            "db",
            "phage_taxon",
            "last_host_taxon",
        ],
        dropna=False,
    )
    .median()
    .sort_index()
)

# Save
df.to_excel(f"{out_folder}/spacers_blast.xlsx")


##################################################
# ---- Blast left phage against OVD and GPD ---- #
##################################################

# Run BLAST
query = left_phage
db = merged_blast_db
out = f"{out_folder}/leftcross__vs__ovd_and_gpd.tsv"  # "/home/usuario/Proyectos/Results/Kelly/Parvimonas/pgap_pipeline/manual/KellyParvimonas_crispr/cross_blastn/leftcross__vs__ovd_and_gpd.tsv"
outfmt = OUTFMT

cmd = [
    "blastn",
    f"-query {query}",
    f"-out {out}",
    f"-db {db}",
    "-perc_identity 50 -max_hsps 5 -max_target_seqs 1000 -num_threads 20",
    f"-outfmt '{outfmt}'",
]
cmd = " ".join(cmd)
print(cmd)
run(cmd, shell=True, check=True)

# Get metadata
merged_ovd_gpd  # comes from previous section

processed_ovd = pd.merge(
    ovd_meta_main_df, ovd_meta_HP_df, on="vOTU ID"
)  # new processed_ovd
processed_ovd = processed_ovd[
    [
        "vOTU ID",
        "Prediction_method",
        "quality__checkv",
        "completeness__checkv",
        "Length (bp)",
    ]
]

# Get results and merge with metadata
res_tmp = pd.read_csv(out, sep="\t", header=None)
res_tmp.columns = [i for i in outfmt.split(" ")[1:]]
merged = pd.merge(merged_ovd_gpd, res_tmp, left_on="sseqid", right_on="sseqid")
merged = pd.merge(merged, processed_ovd, left_on="sseqid", right_on="vOTU ID")
merged = merged.drop_duplicates()

# See top hits
merged = merged.sort_values(["bitscore"], ascending=False)
merged.head(30)
merged.head(30)["last_host_taxon"].unique()


# merged_ovd_gpd[
#     (merged_ovd_gpd["db"] == "oral")
#     & (merged_ovd_gpd["last_host_taxon"].str.contains("Parvimonas"))
# ]


##################################################
# ---------- Pairwise CRISPR blastn ------------ #
##################################################

# Make DB with spacers
cmd = f"makeblastdb -in {merged_spacers} -dbtype nucl"
print(cmd)
run(cmd, shell=True, check=True)

# BLAST spacers against themselves
output_folder = "output"
db = merged_spacers
query = db
out = f"{output_folder}/spacer_pairwise.tsv"
threads = 60
outfmt = "6 qseqid sseqid stitle sstrand pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qseq sseq"

cmd = f"blastn -query {query} \
    -out {out} \
    -db {db} \
    -perc_identity 90 \
    -max_hsps 1 \
    -max_target_seqs 1000 \
    -task blastn-short \
    -num_threads {threads} -outfmt '{outfmt}'\n"

print(cmd)
run(cmd, shell=True, check=True)

# Read results
df = pd.read_csv(out, sep="\t", header=None)
df.columns = [i for i in outfmt.split(" ")[1:]]

# Remove reciprocal hits
df[["qseqid_bc", "qseqid_crispr"]] = df["qseqid"].str.split("__", expand=True)
df[["sseqid_bc", "sseqid_crispr"]] = df["sseqid"].str.split("__", expand=True)
df = df[df["qseqid_bc"] != df["sseqid_bc"]]

# Coverage
df["cov"] = 100 * df["length"] / df["slen"]
df = df[df["cov"] >= 90]

# Sort
df = df.sort_values(by=["sseqid", "cov", "pident"], ascending=False)
df = df[["qseqid", "sseqid", "pident", "cov", "qseq", "sseq"]]

# View
df
