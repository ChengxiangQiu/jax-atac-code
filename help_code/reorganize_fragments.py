

###################################################################
### Reorganize fragments into three levels of cell type annotations

### The fragments from individual samples can be downloaded from GEO
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE325776
### E10.0.fragments.txt.gz
### E10.25.fragments.txt.gz
### ...
### P0.fragments.txt.gz

import gzip
import pandas as pd
from contextlib import ExitStack

work_path = "your_path"
regroup_by = "celltype_L1" ### or celltype_L2, celltype_L3, based on what you need

df = pd.read_csv("https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/df_cell.csv", index_col=0)

sample_cell = df.groupby("day").apply(lambda x: set(x.index)).to_dict()
cell_anno = df[regroup_by].to_dict()
anno_list = list(df[regroup_by].unique())

day_to_num = {
    'E10.25': 1, 'E10.5': 2, 'E11.0': 3, 'E11.75': 4,
    'E12.25': 5, 'E12.5': 6, 'E12.75': 7, 'E13.5': 8,
    'E10.0': 9, 'E10.75': 10, 'E11.25': 11, 'E11.5': 12,
    'E12.0': 13, 'E13.0': 14, 'E13.25': 15, 'E13.75': 16,
    'E14.0': 17, 'E14.25': 18, 'E14.375': 19, 'E14.75': 20,
    'E15.0': 21, 'E15.25': 22, 'E15.5': 23, 'E15.75': 24,
    'E16.0': 25, 'E16.25': 26, 'E16.5': 27, 'E16.75': 28,
    'E17.0': 29, 'E17.25': 30, 'E17.5': 31, 'E18.0': 32,
    'E18.25': 33, 'E18.5': 34, 'E18.75': 35, 'P0': 36
}

with ExitStack() as stack:
    files = {name: stack.enter_context(gzip.open(f'{name}.fragments.txt.gz', 'wt')) for name in anno_list}
    for sample in sample_cell:
        with gzip.open(f"{work_path}/{sample}.fragments.txt.gz", "rt") as file:
            for line in file:
                l = line.rstrip().split('\t')
                cell_id = f"{l[3]}_{day_to_num[sample]}"
                if cell_id in sample_cell[sample]:
                    files[cell_anno[cell_id]].write(line)
