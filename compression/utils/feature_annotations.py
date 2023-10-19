def collect_gene_annotations(anno_fn, genes):
    '''Collect gene annotations from GTF file'''
    if anno_fn is False:
        return None

    featype = 'gene'

    with gzip.open(anno_fn, 'rt') as gtf:
        gene_annos = []
        for line in gtf:
            if f'\t{featype}\t' not in line:
                continue
            fields = line.split('\t')
            if fields[2] != featype:
                continue
            attrs = fields[-1].split(';')

            gene_name = None
            transcript_id = None
            for attr in attrs:
                if 'gene_name' in attr:
                    gene_name = attr.split(' ')[-1][1:-1]
                elif " gene " in attr:
                    gene_name = attr.split(' ')[-1][1:-1]
                elif "ID=gene:" in attr:
                    gene_name = attr.split(';')[0].split(':')[-1]
                # for drosophila, flybase uses "gene_symbol"
                elif 'gene_symbol' in attr:
                    gene_name = attr.split(" ")[-1][1:-1]
                elif 'transcript_id' in attr:
                    transcript_id = attr.split(' ')[-1][1:-1]
                elif 'Name=' in attr:
                    gene_name = attr.split('=')[1]
                    transcript_id = gene_name

            if (gene_name is not None) and (transcript_id is None):
                transcript_id = gene_name

            if (gene_name is None) or (transcript_id is None):
                continue
            gene_annos.append({
                'transcript_id': transcript_id,
                'gene_name': gene_name,
                'chromosome_name': fields[0],
                'start_position': int(fields[3]),
                'end_position': int(fields[4]),
                'strand': 1 if fields[6] == '+' else -1,
                'transcription_start_site': int(fields[3]) if fields[6] == '+' else int(fields[4]),
                })
    gene_annos = pd.DataFrame(gene_annos)

    # NOTE: some species like zebrafish don't really have a transcript id (yet?)
    #assert gene_annos['transcript_id'].value_counts()[0] == 1

    # FIXME: choose the largest transcript or something. For this particular
    # repo it's not that important
    gene_annos = (gene_annos.drop_duplicates('gene_name')
                            .set_index('gene_name', drop=False))

    genes_missing = list(set(genes) - set(gene_annos['gene_name'].values))
    gene_annos_miss = pd.DataFrame([], index=genes_missing)
    gene_annos_miss['transcript_id'] = gene_annos_miss.index
    gene_annos_miss['start_position'] = -1
    gene_annos_miss['end_position'] = -1
    gene_annos_miss['strand'] = 0
    gene_annos_miss['chromosome_name'] = ''
    gene_annos_miss['transcription_start_site'] = -1
    gene_annos = pd.concat([gene_annos, gene_annos_miss])
    gene_annos = gene_annos.loc[genes]
    gene_annos['strand'] = gene_annos['strand'].astype('i2')

    return gene_annos


def collect_feature_annotations(anno_fn, features, measurement_type="gene_expression"):
    if measurement_type == "gene_expression":
        return collect_gene_annotations(anno_fn, features)
    return None



def store_compressed_feature_annotations(
    fn_out,
    feature_annos,
    measurement_type,
    ):
    '''Store compressed feature annotations.'''
    with h5py.File(fn_out, 'a') as h5_data:
        me = h5_data[measurement_type]

        group = me.create_group('feature_annotations')
        group.create_dataset(
                'feature_name',
                data=feature_annos.index.values.astype('S'))
        if 'transcription_start_site' in feature_annos.columns:
            group.create_dataset(
                    'transcription_start_site',
                    data=feature_annos['transcription_start_site'].values, dtype='i8')
        group.create_dataset(
                'chromosome_name',
                data=feature_annos['chromosome_name'].astype('S'))
        group.create_dataset(
                'start_position',
                data=feature_annos['start_position'].values, dtype='i8')
        group.create_dataset(
                'end_position',
                data=feature_annos['end_position'].values, dtype='i8')
        group.create_dataset(
                'strand', data=feature_annos['strand'].values, dtype='i2')
