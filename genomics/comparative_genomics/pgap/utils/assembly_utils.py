from Bio import SeqIO
import re

def rename_contigs(input_assembly, output_assembly, wc):
    """
    Rename contigs in assembly

    Rules:
    - Contig ID is just a number: replace by {barcode}.{id}
    - Contig ID fits contig_\d+: replace by {barcode}.{id}, {id} being the number in contig_\d+
    - Contig ID contains barcode: do nothing
    - Else: do nothing

    Input:
    - Assembly path in fasta format
    Output:
    - Assembly path in fasta format
    """
    with open(output_assembly, "wt") as handle:
        for record in SeqIO.parse(input_assembly, "fasta"):
            if re.findall(r'^\d.*',record.id):
                # If ID starts with numbers (unicycler outputs contigs like this: ">1" )
                # CHANGE: barcode + '.' + id
                record.id = f'{wc.barcode}.{record.id}'

            elif re.findall(r'^contig_\d.*',record.id):
                # Flye outputs contigs names like "contig_N"
                # CHANGE: barcode + the contig number
                record.id = record.id.split('_')[1]
                record.id = f'{wc.barcode}.{record.id}'

            elif wc.barcode in record.id:
                # If the barcode is already inside dont change
                pass

            else:
                # Else dont touch the contig id because it probably has an unique id
                pass
            
            SeqIO.write(record, handle, "fasta")
