"""Per run: top-level features in a ProteomicsLFQ -feat_dir checkpoint, with and without peptide IDs."""
import sys, glob, os, collections
import pyarrow.parquet as pq
for d in sorted(glob.glob(os.path.join(sys.argv[1], '*.featureParquet'))):
    f = pq.read_table(os.path.join(d, 'features.parquet'), columns=['unique_id', 'depth', 'charge']).to_pandas()
    top = f[f.depth == 0]
    p = pq.read_table(os.path.join(d, 'psms.parquet'), columns=['feature_unique_id']).to_pandas()
    with_id = top.unique_id.isin(set(p.feature_unique_id))
    print(os.path.basename(d).split('_fg')[0], len(top), int(with_id.sum()), int((~with_id).sum()),
          dict(collections.Counter(top[~with_id].charge)), sep='\t')
