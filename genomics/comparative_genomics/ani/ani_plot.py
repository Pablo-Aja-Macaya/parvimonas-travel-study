import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import pandas as pd
from io import StringIO
import os


def process_orthoaniu_file(ani_list_file: str) -> pd.DataFrame:
    """
    Takes in a file which is the result of running OrthoAniU
    in "list" format, and separates the result into matrix and list format
    """

    with open(ani_list_file, "rt") as input_ani:
        lines = input_ani.readlines()

    list_start = lines.index("# Pair number mappings\n")
    matrix_start = lines.index(
        "pair#\torthoANI_value\tavg_aligned_length\tquery_coverage\tsubject_coverage\tquery_length\tsubject_length\n"
    )

    ani_list = lines[list_start + 1 : matrix_start - 2]
    matrix_list = lines[matrix_start:]

    ani_list = StringIO(
        "\n".join(ani_list).replace("#", "").replace("|", "\t").replace(" ", "")
    )
    matrix_list = StringIO("\n".join(matrix_list).replace("#", "").replace(" ", ""))

    # ---- Process ANI results ----
    df_pairs = pd.read_csv(ani_list, sep="\t")
    df_matrix = pd.read_csv(matrix_list, sep="\t")

    df_pairs.columns = ["sample1", "sample2"]

    # Concat matrix and names
    df_c = pd.concat([df_pairs.reset_index(drop=True), df_matrix], axis=1)

    # Order
    df_c = df_c.sort_values("orthoANI_value", ascending=False)

    # Samples
    unique_samples = pd.unique(df_c[["sample1", "sample2"]].values.ravel())
    rows = []
    for sample in unique_samples:
        # Load
        direct = df_c[df_c["sample1"] == sample][["sample2", "orthoANI_value"]]
        indirect = df_c[df_c["sample2"] == sample][["sample1", "orthoANI_value"]]

        # Rename
        direct.columns = ["sample", "orthoANI_value"]
        indirect.columns = ["sample", "orthoANI_value"]

        # Concat and order
        df = pd.concat([indirect, direct], axis=0)
        df = df.append({"sample": sample, "orthoANI_value": 100}, ignore_index=True)
        df = df.sort_values("sample", ascending=False)

        names = []
        row = []
        for s, value in df.values:
            names.append(s)
            row.append(value)
        rows.append(row)

    matrix_df = pd.DataFrame(rows, index=unique_samples, columns=names)
    matrix_df = matrix_df.sort_index(ascending=False)

    return df_c, matrix_df


# ---- Input ----
# OrthoANIu result path using "list" format (Instal from https://www.ezbiocloud.net/tools/orthoaniu)
# Command: java -jar $ortho_ani_u_path -u $usearch_path -fd $input_folder -o $output_res -fmt list -n 12
ani_result = "ani_list_output.txt"

# Names to change row and column wise
corr_dicc = {
    "GCF_016552165.1": {"strain": "S3374"},
    "GCA_000154405.1": {"strain": "ATCC 33270"},
    "GCA_000493795.1": {"strain": "A293"},
    "GCA_000800295.1": {"strain": "KCOM 1535"},
    "GCA_003454775.1": {"strain": "KCOM 1037"},
    "GCA_003938725.1": {"strain": "FDAARGOS_569"},
    "GCA_008633155.1": {"strain": "13-07-26"},
    "GCA_900637905.1": {"strain": "NCTC11808"},
    "GCA_902373675.1": {"strain": "MGYG-HGUT-01301"},
    "GCA_023147965.1": {"strain": "EYE_66"},
    "GCA_023148065.1": {"strain": "EYE_30"},
    "GCA_023148105.1": {"strain": "EYE_29"},
    "GCA_023148225.1": {"strain": "EYE_25"},
    "PM102KC_G_1": {"strain": "PM102KC-G-1"},
    "PM114KC_G_1": {"strain": "PM114KC-G-1"},
    "PM79KC_AC_1": {"strain": "PM79KC-AC-1"},
    "PM79KC_AC_2": {"strain": "PM79KC-AC-2"},
    "PM79KC_AC_3": {"strain": "PM79KC-AC-3"},
    "PM79KC_AC_4": {"strain": "PM79KC-AC-4"},
    "PM79KC_G_1": {"strain": "PM79KC-G-1"},
    "PM89KC_AC_1": {"strain": "PM89KC-AC-1"},
    "PM89KC_G_1": {"strain": "PM89KC-G-1"},
    "PM89KC_G_2": {"strain": "PM89KC-G-2"},
    "PM94KC_G_1": {"strain": "PM94KC-G-1"},
}


# ---- Process and plot heatmap ----
list_df, matrix_df = process_orthoaniu_file(ani_result)


# Format matrix dataframe
matrix_df = matrix_df.rename(
    columns={f"{k}.fna": v["strain"] for k, v in corr_dicc.items()}
)
matrix_df = matrix_df.rename({f"{k}.fna": v["strain"] for k, v in corr_dicc.items()})
matrix_df.index.names = [None]

# Plot
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"
g = sns.clustermap(
    matrix_df,
    cmap="RdYlGn",
    vmin=96,
    vmax=100,
    metric="euclidean",
    annot=True,
    fmt=".2f",
    annot_kws={"fontsize": 8},
    cbar_kws={"label": "OrthoANI"},
)
g.fig.subplots_adjust(right=0.7)
g.ax_cbar.set_position((0.8, 0.2, 0.03, 0.4))
g.fig.set_figwidth(20)
g.fig.set_figheight(11.7)

# Save figures
plt.savefig(
    "ani.png",
    dpi=300,
)
plt.savefig(
    "ani.tiff",
    dpi=300,
)
# plt.show()
