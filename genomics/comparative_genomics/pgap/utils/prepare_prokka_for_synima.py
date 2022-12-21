import glob
import shutil
import os
import re

def prepare_prokka_for_synima(barcode, input_folder, output_folder):
    """
    Prepare Prokka's .faa .ffn .gff and .fna for Synima (Synteny analysis)
    Includes code from https://github.com/Ivan-Pchelin
    """

    # ---- Settings ----
    # barcode = '4i3-ccr89tumor_circlator'
    # input_folder= f'output/pruebas/prokka_transformation/org/{barcode}'
    # output_folder = f'output/pruebas/prokka_transformation/mod/{barcode}'

    compliant_barcode = barcode.replace('-','_')
    
    # ---- Create output folder ----
    os.makedirs(output_folder, exist_ok=True)

    # ---- Change suffixes ----
    prokka_suffix_corr = {
        'faa' : 'annotation.pep',
        'ffn' : 'annotation.cds',
        'gff' : 'annotation.gff3',
        'fna' : 'genome.fa'
    }
    mod_files = []
    for org_suffix,new_suffix in prokka_suffix_corr.items():
        f = shutil.copy2(f'{input_folder}/{barcode}.{org_suffix}', f'{output_folder}/{compliant_barcode}.{new_suffix}')
        mod_files.append(f)

    # ---- Modify files to fit Synima ----
    # (https://github.com/Ivan-Pchelin/scripts-for-synteny-visualization/blob/main/4_prepare_for_synima.py)

    sequencenames = set()
    files = glob.glob(f'{output_folder}/*.*', recursive=True)
    for f in files:
        if f.endswith('cds') or f.endswith('pep'):
            with open (f) as inf:
                lines = inf.readlines()
            with open (f, 'w') as ouf:
                for i in lines:
                    if '>' not in i:
                        ouf.write(i)
                    else:
                        if i[0:5] == '>cds-':
                            ouf.write(i)
                        else:
                            ouf.write('>cds-'+i[1:])
        if f.endswith('gff3'):
            with open (f) as inf:
                lines = inf.readlines()
            with open (f, 'w') as ouf:
                filename = os.path.basename(f)
                sequencename = re.split(r'\.',filename)[-3]
                sequencenames.add(sequencename)
                for i in lines:
                    if len(re.findall(r'\w+', i)) > 1:
                        if '#' not in i and '>' not in i:
                            if re.findall('[^ATGC]', i.strip()):
                                stuff = re.split(r'\t', i)
                                stuff[1] = sequencename
                                stuff[2] = 'gene'
                            for i in stuff:
                                ouf.write(i)
                                if '\n' not in i:
                                    ouf.write('\t')

    # Create a global set with entries from PEP files
    pep_files = glob.glob(f'{output_folder}/*.pep', recursive=True)
    cds_names = set()
    for i in pep_files:
        with open (i) as inf:
            lines = inf.readlines()
        for l in lines:
            if l.startswith('>cds-'):
                cds_names.add(re.split(r' ', l)[0][5:])

    # Delete from CDS files all entries missing from PEP files
    for f in files:
        if '.py' not in f and 'Repo_spec' not in f:
            if f.endswith('cds'):
                with open (f) as inf:
                    allseqs = []
                    rawlines = inf.readlines()
                with open (f, 'w') as ouf:
                    lines = []
                    applicant = ''
                    letts = ''
                    d = 0
                    for g in rawlines:
                        if re.findall(r'\w', g):
                            if '>' in g:
                                if d != 0:
                                    allseqs.append(applicant.strip())
                                    allseqs.append(letts)
                                    letts = ''
                                applicant = g
                                d += 1
                            else:
                                letts += g.strip()
                    allseqs.append(applicant.strip())
                    allseqs.append(letts)

                    runner = 0                    
                    while runner < len(allseqs):
                        if allseqs[runner].startswith('>'):
                            if re.findall(r'\S+', allseqs[runner])[0][5:] in cds_names:
                                ouf.write(allseqs[runner])
                                ouf.write('\n')
                                ouf.write(allseqs[runner+1])
                                ouf.write('\n')
                        runner += 2
    print('\n CDS files are done')
    
    # Delete from GFF3 files all entries missing from PEP files
    gff_files = glob.glob(f'{output_folder}/*.gff3', recursive=True)
    print(gff_files)
    for gff in gff_files:
        with open (gff) as inf:
            rawlines = inf.readlines()
        with open (gff, 'w') as ouf:
            for raw in rawlines:
                if len(re.findall(r'\w', raw)) > 2:
                    for cds in cds_names:
                        
                        if cds in re.split(r';', raw)[0]:
                            outentry = re.split(r'\t', raw)
                            
                            if outentry[8][0:7] != 'ID=cds-':
                                outentry[8] = 'ID=cds-' + outentry[8][3:]
                                
                            a = 0
                            while a < len(outentry):
                                ouf.write(outentry[a])
                                if a < len(outentry) - 1:
                                    ouf.write('\t')
                                a += 1
    print('\n GFF3 files are done')
    
    return mod_files