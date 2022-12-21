
# Rename Kelly's tree

f = 'output/NanoporeBacteria/Kelly/pgap_pipeline/pipeline/comparative_analysis/03_pirate/output_folder/binary_presence_absence.nwk'
with open(f, 'rt') as handle:
    tree = handle.readlines()[0]

# c("cluster_order", "ATCC 33270", "A293", "KCOM 1535", 
#                              "KCOM 1037", "FDAARGOS_569", "13-07-26", "EYE_66", 
#                              "EYE_30", "EYE_29", "EYE_25", "NCTC11808", 
#                              "MGYG-HGUT-01301", "S3374", "PM102KC-G-1", "PM114KC-G-1", 
#                              "PM79KC-AC-1", "PM79KC-AC-2", "PM79KC-AC-3", "PM79KC-AC-4", "PM79KC-G-1", 
#                              "PM89KC-AC-1", "PM89KC-G-1", "PM89KC-G-2", "PM94KC-G-1")
corr_dicc = {
    "GCF_016552165_1":{'cepa_name':'S3374'},
    "GCA_000154405_1":{'cepa_name':'ATCC 33270'},
    "GCA_000493795_1":{'cepa_name':'A293'},
    "GCA_000800295_1":{'cepa_name':'KCOM 1535'},
    "GCA_003454775_1":{'cepa_name':'KCOM 1037'},
    "GCA_003938725_1":{'cepa_name':'FDAARGOS_569'},
    "GCA_008633155_1":{'cepa_name':'13-07-26'},
    "GCA_900637905_1":{'cepa_name':'NCTC11808'},
    "GCA_902373675_1":{'cepa_name':'MGYG-HGUT-01301'},
    "GCA_023147965_1":{'cepa_name':'EYE_66'},
    "GCA_023148065_1":{'cepa_name':'EYE_30'},
    "GCA_023148105_1":{'cepa_name':'EYE_29'},
    "GCA_023148225_1":{'cepa_name':'EYE_25'},
    "PM102KC_G_1":{'cepa_name':'PM102KC-G-1'},
    "PM114KC_G_1":{'cepa_name':'PM114KC-G-1'},
    "PM79KC_AC_1":{'cepa_name':'PM79KC-AC-1'},
    "PM79KC_AC_2":{'cepa_name':'PM79KC-AC-2'},
    "PM79KC_AC_3":{'cepa_name':'PM79KC-AC-3'},
    "PM79KC_AC_4":{'cepa_name':'PM79KC-AC-4'},
    "PM79KC_G_1":{'cepa_name':'PM79KC-G-1'},
    "PM89KC_AC_1":{'cepa_name':'PM89KC-AC-1'},
    "PM89KC_G_1":{'cepa_name':'PM89KC-G-1'},
    "PM89KC_G_2":{'cepa_name':'PM89KC-G-2'},
    "PM94KC_G_1":{'cepa_name':'PM94KC-G-1'},
}


for k,v in corr_dicc.items():
    # k = k.replace('-','_')
    tree = tree.replace(k,v['cepa_name'])

with open('output/NanoporeBacteria/Kelly/pgap_pipeline/pipeline/comparative_analysis/03_pirate/output_folder/binary_presence_absence.renamed.nwk','wt') as handle:
    handle.write(tree)

