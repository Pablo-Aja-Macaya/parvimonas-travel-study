corr = {
    '1b8-ccr79tumor':	'barcode01',
    '2d4-ccr79normal':	'barcode02',
    '1c4-ccr79tumor':	'barcode03',
    '2g7-ccr79oral':	'barcode04',
    '1d8-ccr79tumor':	'barcode05',
    '1b5-ccr79tumor':	'barcode06',
    '5b8-ccr94oral':	'barcode07',
    '6a5-ccr89oral':	'barcode08',
    '6a7-ccr89oral':	'barcode09',
    '6e2-ccr102oral':	'barcode10',
    '4i3-ccr89tumor':	'barcode11',
    '7f8-ccr114tumor':	'barcode12',    
}

# To illumina
input_folder = 'output/NanoporeBacteria/Kelly/pmicra_merged__flye/reads_qc/07_clean_reads'
output_folder = 'output/NanoporeBacteria/Kelly/pmicra_merged__unicycler_hybrid/00_long_clean_reads'
for illumina_id, barcode in corr.items():
    print(f'mkdir -p {output_folder}/{illumina_id}')
    print(f'cp {input_folder}/{barcode}/{barcode}.fastq.gz {output_folder}/{illumina_id}/{illumina_id}.fastq.gz' )


# From illumina
output_folder = 'output/NanoporeBacteria/Kelly/pmicra040322__flye__readsplitting/reads_qc/07_clean_reads'
input_folder = 'output/NanoporeBacteria/Kelly/pmicra040322__unicycler_hybrid__readsplitting/00_long_clean_reads'
for illumina_id, barcode in corr.items():
    print(f'cp {input_folder}/{illumina_id}/{illumina_id}.fastq.gz {output_folder}/{barcode}/{barcode}.fastq.gz' )
