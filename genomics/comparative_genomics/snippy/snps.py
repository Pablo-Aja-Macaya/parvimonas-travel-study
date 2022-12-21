import itertools
import pandas as pd
from subprocess import run

# Snippy version used: snippy 4.6.0

# Input
input_folder = "/home/usuario/Proyectos/Innova/output/NanoporeBacteria/Kelly/pgap_pipeline/pipeline/03_converted_to_prokka"
main_output_folder = "/home/usuario/Proyectos/Innova/output/NanoporeBacteria/Kelly/pgap_pipeline/manual/snippy"

sample_dict = {
    "PM89KC_G_1": {
        "assembly": f"{input_folder}/PM89KC_G_1/PM89KC_G_1.fna",
        "annotation": f"{input_folder}/PM89KC_G_1/PM89KC_G_1.gbk",
    },
    "PM89KC_G_2": {
        "assembly": f"{input_folder}/PM89KC_G_2/PM89KC_G_2.fna",
        "annotation": f"{input_folder}/PM89KC_G_2/PM89KC_G_2.gbk",
    },
    "PM89KC_AC_1": {
        "assembly": f"{input_folder}/PM89KC_AC_1/PM89KC_AC_1.fna",
        "annotation": f"{input_folder}/PM89KC_AC_1/PM89KC_AC_1.gbk",
    },
}

# ---- Run snippy ----
# Pick every sample as reference and compare to the rest with snippy
threads = 60
snippy_results = {}
for ref_key, query_key in itertools.combinations(sample_dict, 2):
    print(f"\n# {ref_key} vs {query_key}")
    reference = sample_dict[ref_key]["annotation"]
    query = sample_dict[query_key]["assembly"]
    output_folder = f"{main_output_folder}/ref_{ref_key}__vs__{query_key}"

    cmd = f"snippy --ctgs {query} --ref {reference} --outdir {output_folder} --cpus {threads} --cleanup"
    _ = run(cmd, shell=True, check=True)

    if snippy_results.get(ref_key):
        snippy_results[ref_key] += [
            {"query": query_key, "results_folder": output_folder}
        ]
    else:
        snippy_results[ref_key] = [
            {"query": query_key, "results_folder": output_folder}
        ]

# ---- Get stats ----
for ref, queries in snippy_results.items():
    for query in queries:
        snps_file = f"{query['results_folder']}/snps.tab"

        # No filtering
        snps_df = pd.read_csv(snps_file, sep="\t")
        full_snps_counts = snps_df["TYPE"].value_counts().to_dict()

        # Filtering to remove non synonymous
        non_synonymous_snps_df = snps_df[
            ~(snps_df["EFFECT"].str.contains("synonymous", na=False))
        ]
        ns_snps_counts = non_synonymous_snps_df["TYPE"].value_counts().to_dict()

        # Print
        with open(f"{query['results_folder']}/res.txt", "wt") as handle:
            handle.write(f'# ref {ref} vs {query["query"]}\n')
            handle.write(f"All mutations: {full_snps_counts}\n")
            handle.write(f"NS  mutations: {ns_snps_counts}\n")

        # Output non synonymous mutations per gene
        grouped_muts = pd.DataFrame(
            non_synonymous_snps_df.groupby(
                by=["LOCUS_TAG", "GENE", "PRODUCT"], dropna=False
            ).size()
        )
        grouped_muts.columns = ["COUNT"]
        grouped_muts = grouped_muts.reset_index()  # reset to table mode
        grouped_muts = grouped_muts[
            ~grouped_muts["LOCUS_TAG"].isna()
        ]  # remove the row without LOCUS_TAG
        grouped_muts = grouped_muts.sort_values("COUNT", ascending=False)
        grouped_muts.to_csv(
            f"{query['results_folder']}/non_synonymous_snps_counts.tsv",
            sep="\t",
            index=None,
        )
