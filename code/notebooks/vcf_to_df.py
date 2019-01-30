import glob
import pandas as pd
from pathlib import Path
from pysam import VariantFile


def vcf_to_dataframe(vcf_file_path):


    path = Path(filename)
    sample = path.stem

    vcf_in = VariantFile(filename)
    vcf_out = VariantFile('{}/{}.filter.vcf'.format(str(path.parent),str(path.stem)),
                          'w', header=vcf_in.header)

    lst_ = []

    for rec in vcf_in:
        rec = rec
        entry = [rec.contig,
                rec.start,
                rec.stop,
                rec.ref,
                rec.id,
                rec.alts[0],
                ",".join(rec.filter.keys()),
                rec.qual,
                rec.info.get('QUAL'),
                rec.info.get('AF'),
                rec.info.get('DP'),
                rec.info.get('VD'),
                rec.info.get('TYPE'),
                rec.info.get('PMEAN'),
                rec.info.get('NM'),
                CAF_to_MAF(rec.info.get('CAF'))]

        chr_,s,e,ref,id_,alt,filter_,qual,qual_,af,dp,vd,type_,pmean,nm,GMAF = entry

        if ((str(filter_) == filters['filter']) & (str(type_).isin(filters['type'])):
            if (float(af) > filters['af']) & (pmean > filters['pmean']):
                if not GMAF or float(GMAF) < filters['gmaf']:
                    if int(dp) > 10:
                        lst_.append(entry)
                        vcf_out.write(rec)

    df = pd.DataFrame(np.row_stack(lst_),columns=['chr',
                                                  'start',
                                                  'stop',
                                                  'ref',
                                                  'id',
                                                  'alt',
                                                  'filter',
                                                  'qual',
                                                  'vd_qual',
                                                  'af',
                                                  'dp',
                                                  'vd',
                                                  'type',
                                                  'pmean',
                                                  'nm',
                                                  'GMAF'])

    df['af'] = df['af'].astype('float')
    df['pmean'] = df['pmean'].astype('float')
    df['dp'] = [int(dp) if dp else None for dp in df['dp']]
    df['GMAF'] = [float(MAF) if MAF else 0 for MAF in df['GMAF']]

    if sample.split('.')[0] == 'FRFZ':
        df['chr'] = ['chr' + c for c in df['chr']]

    df.to_pickle(str(path.parent) + '/pkls/' + str(path.stem) + '.pkl')
