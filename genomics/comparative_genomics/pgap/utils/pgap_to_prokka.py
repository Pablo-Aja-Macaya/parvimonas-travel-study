from Bio import SeqIO

# ---- GFF transformation ----
def pgap_gff_to_prokka_gff(barcode, input_assembly, input_assembly_gff, output_assembly_gff):
    """
    Transform PGAP's GFF format to Prokka's GFF format
    Main difference is PGAP's GFF does not contain the contig sequences
    """
    compliant_barcode = barcode.replace('-','_')

    # Process GFF
    with open(input_assembly_gff, 'rt') as input_handle:
        sequence_regions = []
        lines_to_write = []
        for line in input_handle:
            line = line.strip()
            line = line.replace('>lcl|', '>')
            line = line.replace(barcode, compliant_barcode)
            if '#!processor NCBI annotwriter' in line or '#!gff-spec-version 1.21' in line or '##species' in line or '##gff-version 3' in line:
                # Dont keep it
                pass
            elif '##sequence-region' in line:
                #Extract region (##sequence-region start_n end_n)
                sequence_regions.append(line)

            else:
                # Extrac main components of GFF
                gff_atrs = {}
                gff_atrs['seqname'], gff_atrs['source'], gff_atrs['feature'], gff_atrs['start'], gff_atrs['end'], gff_atrs['score'], gff_atrs['strand'], gff_atrs['frame'], gff_atrs['attribute']  = line.split('\t')
                
                # Create dict for atributes
                tmp_atrs = {a.split('=')[0]:a.split('=')[1] for a in gff_atrs['attribute'].split(';')}
                gff_atrs['attribute'] = tmp_atrs

                # If the atribute is a certain type
                if gff_atrs['feature'] != 'gene':
                    # Replace ids
                    gff_atrs['attribute']['ID'] = gff_atrs['attribute']['ID'].replace('cds-','')

                    if gff_atrs['attribute'].get('Parent'):
                        gff_atrs['attribute']['Parent'] = gff_atrs['attribute']['Parent'].replace('gene-','')

                    # Return dict to string in GFF format
                    gff_atrs['attribute'] = ';'.join([f'{k}={v}' for k,v in gff_atrs['attribute'].items()])

                    lines_to_write.append('\t'.join(gff_atrs.values()))

                else:
                    pass

    # Write GFF data to file
    with open(output_assembly_gff,'wt') as output_handle:
        output_handle.write('##gff-version 3\n')   
        for i in sequence_regions:
            output_handle.write(f'{i}\n')   
        for i in lines_to_write:
            output_handle.write(f'{i}\n')

    # Write contig data to file
    with open(output_assembly_gff, "at") as output_handle:
        output_handle.write(f'##FASTA\n')
        for record in SeqIO.parse(input_assembly, "fasta"):
            record.id = record.id.replace('lcl|', '').replace('-','_')
            record.description = ''
            record.name = ''
            SeqIO.write(record, output_handle, "fasta")
    
    return output_assembly_gff

# ---- File with all genes in nucleotide form ----
def extract_ncl_seqs_from_gbk(gbk_input, ncl_genes_output):
    """
    Extract nucleotide sequences in genbank file annotated with PGAP to 
    an output file, creating Prokka's .ffn file
    """
    # Reset output file if it exists
    with open(ncl_genes_output, 'wt') as handle:
        pass

    # For every contig
    for contig in SeqIO.parse(gbk_input, "genbank"):
        # Get contig features and for every feature extract relevant information
        source = contig.annotations['source']
        if contig.features:
            for feature in contig.features:
                    # Type of feature 1
                    if feature.type in ['CDS','rRNA','tRNA','tmRNA']:
                        qualifiers = feature.qualifiers
                        ncl_sequence = feature.extract(contig.seq)
                        locus_tag = qualifiers['locus_tag'][0]
                        product = qualifiers['product'][0]
                        
                        with open(ncl_genes_output, 'at') as handle:
                            handle.write(f'>{locus_tag} {product} [{source}]\n{ncl_sequence}\n')

                    # Type of feature 2
                    elif feature.type in ['misc_binding','misc_feature','regulatory']:
                        qualifiers = feature.qualifiers
                        ncl_sequence = feature.extract(contig.seq)
                        locus_tag = 'no_locus_tag'
                        product = qualifiers['note'][0].split(';')[0]

                        # print(f'>{locus_tag} {product} [{source}]\n{ncl_sequence}')

                    # Type of feature 3
                    elif feature.type in ['repeat_region']:
                        qualifiers = feature.qualifiers
                        ncl_sequence = feature.extract(contig.seq)
                        locus_tag = 'no_locus_tag'
                        product = qualifiers['rpt_family'][0]
                        repeat = qualifiers['rpt_unit_seq']

                        # print(f'>{locus_tag} {product} [{source}]\n{ncl_sequence}')

                    elif feature.type in ['ncRNA']:
                        qualifiers = feature.qualifiers
                        ncl_sequence = feature.extract(contig.seq)
                        locus_tag = qualifiers['locus_tag'][0]
                        product = qualifiers['product'][0]                  

                        # print(f'>{locus_tag} {product} [{source}]\n{ncl_sequence}')
    
    return ncl_genes_output

# ---- Genome and proteins ----
def clean_pgap_fasta(input_file, output_file, version=None):
    """
    Clean ids of FASTA files (.faa and .ffn) to remove certain PGAP stuff
    """
    with open(output_file, "xt") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.description = record.description.replace(f"{record.id} ", '')
            record.name = record.name.replace(f"{record.id}", '')
            record.id = record.id.replace('lcl|', '').replace('gnl|','').replace('extdb|','').replace('-','_')
            record.name = ''
            if version=='fna':
                record.description = ''
                
            SeqIO.write(record, output_handle, "fasta")

