import pandas as pd
import numpy as np
import functools
import pathlib


def new_or_missing_pirate(
    gene_families_file: str, output_folder: str, query_groups: dict
):
    """
    Process PIRATE's gene family file, comparing each sample to samples in the query_group dictionary
    Example of query_group:
        query_groups = { # Names must not contain dots because PIRATE replaces them by underscore
            'group_1':['PM89KC_G_1','PM89KC_G_2'],
            'group_2':['PM79KC_AC_1','PM79KC_AC_2','PM79KC_AC_3','PM79KC_AC_4'],
        }
    """
    try:
        with open(f"{output_folder}/groups.txt", "wt") as handle:
            for k, v in query_groups.items():
                handle.write(f"{k}: {v}\n")

        # Read data
        org_gene_families_df = pd.read_csv(gene_families_file, sep="\t")

        # Select sample cols
        pirate_cols = [
            "allele_name",
            "gene_family",
            "consensus_gene_name",
            "consensus_product",
            "threshold",
            "alleles_at_maximum_threshold",
            "number_genomes",
            "average_dose",
            "min_dose",
            "max_dose",
            "genomes_containing_fissions",
            "genomes_containing_duplications",
            "number_fission_loci",
            "number_duplicated_loci",
            "no_loci",
            "products",
            "gene_names",
            "min_length(bp)",
            "max_length(bp)",
            "average_length(bp)",
            "cluster",
            "cluster_order",
        ]

        all_sample_cols = [
            c for c in org_gene_families_df.columns if c not in pirate_cols
        ]

        print("Comparing samples in PIRATE result...")
        for main_target in all_sample_cols:
            # print()
            # print('~'*70)
            # print(main_target.center(70))
            # print('~'*70)
            for query_group_name, query_group_elements in query_groups.items():

                query_group_elements = [
                    i for i in query_group_elements if i != main_target
                ]
                if not query_group_elements:
                    continue

                # print(f'\nVersus: {query_group_name} ({query_group_elements})')

                sample_cols = [main_target] + query_group_elements

                families_target_cols = pirate_cols + sample_cols

                # Filter columns
                gene_families_df = org_gene_families_df[families_target_cols]

                # ---- Gene families ----
                other_components = [i for i in sample_cols if i != main_target]

                def disjunction(*conditions):
                    """
                    Merge multiple conditions with OR
                    """
                    return functools.reduce(np.logical_or, conditions)

                # Get gene families not in target
                condition = [
                    ~gene_families_df[sample].isna() for sample in other_components
                ]
                fams_not_in_target = gene_families_df[
                    (gene_families_df[main_target].isna()) & disjunction(*condition)
                ]
                # print(f'Families NIT: {len(fams_not_in_target)}')
                # print(fams_not_in_target)

                # Get gene families only in target
                condition = [
                    gene_families_df[sample].isna() for sample in other_components
                ]
                fams_only_in_target = gene_families_df[
                    (~gene_families_df[main_target].isna()) & disjunction(*condition)
                ]
                # print(f'Families OIT: {len(fams_only_in_target)}')
                # print(fams_only_in_target)

                # ---- Save ----
                tmp_folder = f"{output_folder}/{main_target}"
                pathlib.Path(tmp_folder).mkdir(parents=True, exist_ok=True)
                fams_not_in_target.to_csv(
                    f"{tmp_folder}/{query_group_name}__not_in_{main_target}.tsv",
                    sep="\t",
                    index=None,
                )
                fams_only_in_target.to_csv(
                    f"{tmp_folder}/{query_group_name}__only_in_{main_target}.tsv",
                    sep="\t",
                    index=None,
                )

    except Exception as e:
        print(e)
        raise Exception(e)
