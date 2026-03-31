### read mtx to scanpy object

import anndata
import scanpy as sc
import pandas as pd
import sys, os

work_path = ""

sample_id = str(int(sys.argv[1]))   ### 1, ..., 36

atac = sc.read(os.path.join(work_path, "gene_count_%s.mtx"%sample_id))
atac.obs = pd.read_csv(os.path.join(work_path, "df_cell_%s.csv"%sample_id), index_col = 0)
atac.var = pd.read_csv(os.path.join(work_path, "df_gene_%s.csv"%sample_id), index_col = 0)

atac.write(os.path.join(work_path, "atac_%s.h5ad"%sample_id), compression="gzip")

import anndata
import scanpy as sc
import pandas as pd
import sys, os

### read atac dataset

adata_1 = sc.read_h5ad(os.path.join(work_path, "atac_1.h5ad"))
adata_2 = sc.read_h5ad(os.path.join(work_path, "atac_2.h5ad"))
adata_3 = sc.read_h5ad(os.path.join(work_path, "atac_3.h5ad"))
adata_4 = sc.read_h5ad(os.path.join(work_path, "atac_4.h5ad"))
adata_5 = sc.read_h5ad(os.path.join(work_path, "atac_5.h5ad"))

adata_6 = sc.read_h5ad(os.path.join(work_path, "atac_6.h5ad"))
adata_7 = sc.read_h5ad(os.path.join(work_path, "atac_7.h5ad"))
adata_8 = sc.read_h5ad(os.path.join(work_path, "atac_8.h5ad"))
adata_9 = sc.read_h5ad(os.path.join(work_path, "atac_9.h5ad"))
adata_10 = sc.read_h5ad(os.path.join(work_path, "atac_10.h5ad"))

adata_11 = sc.read_h5ad(os.path.join(work_path, "atac_11.h5ad"))
adata_12 = sc.read_h5ad(os.path.join(work_path, "atac_12.h5ad"))
adata_13 = sc.read_h5ad(os.path.join(work_path, "atac_13.h5ad"))
adata_14 = sc.read_h5ad(os.path.join(work_path, "atac_14.h5ad"))
adata_15 = sc.read_h5ad(os.path.join(work_path, "atac_15.h5ad"))

adata_16 = sc.read_h5ad(os.path.join(work_path, "atac_16.h5ad"))
adata_17 = sc.read_h5ad(os.path.join(work_path, "atac_17.h5ad"))
adata_18 = sc.read_h5ad(os.path.join(work_path, "atac_18.h5ad"))
adata_19 = sc.read_h5ad(os.path.join(work_path, "atac_19.h5ad"))
adata_20 = sc.read_h5ad(os.path.join(work_path, "atac_20.h5ad"))

adata_21 = sc.read_h5ad(os.path.join(work_path, "atac_21.h5ad"))
adata_22 = sc.read_h5ad(os.path.join(work_path, "atac_22.h5ad"))
adata_23 = sc.read_h5ad(os.path.join(work_path, "atac_23.h5ad"))
adata_24 = sc.read_h5ad(os.path.join(work_path, "atac_24.h5ad"))
adata_25 = sc.read_h5ad(os.path.join(work_path, "atac_25.h5ad"))

adata_26 = sc.read_h5ad(os.path.join(work_path, "atac_26.h5ad"))
adata_27 = sc.read_h5ad(os.path.join(work_path, "atac_27.h5ad"))
adata_28 = sc.read_h5ad(os.path.join(work_path, "atac_28.h5ad"))
adata_29 = sc.read_h5ad(os.path.join(work_path, "atac_29.h5ad"))
adata_30 = sc.read_h5ad(os.path.join(work_path, "atac_30.h5ad"))

adata_31 = sc.read_h5ad(os.path.join(work_path, "atac_31.h5ad"))
adata_32 = sc.read_h5ad(os.path.join(work_path, "atac_32.h5ad"))
adata_33 = sc.read_h5ad(os.path.join(work_path, "atac_33.h5ad"))
adata_34 = sc.read_h5ad(os.path.join(work_path, "atac_34.h5ad"))
adata_35 = sc.read_h5ad(os.path.join(work_path, "atac_35.h5ad"))

adata_36 = sc.read_h5ad(os.path.join(work_path, "atac_36.h5ad"))

adata = adata_1.concatenate(adata_2, adata_3, adata_4, adata_5, 
                            adata_6, adata_7, adata_8, adata_9, adata_10, 
                            adata_11, adata_12, adata_13, adata_14, adata_15,
                            adata_16, adata_17, adata_18, adata_19, adata_20, 
                            adata_21, adata_22, adata_23, adata_24, adata_25,
                            adata_26, adata_27, adata_28, adata_29, adata_30, 
                            adata_31, adata_32, adata_33, adata_34, adata_35,
                            adata_36)

adata.write(os.path.join(work_path, "atac_all.h5ad"), compression="gzip")
adata.obs.to_csv(os.path.join(work_path, 'atac_all.obs.csv'))
adata.var.to_csv(os.path.join(work_path, 'atac_all.var.csv'))


### the h5ad object, "atac_all.h5ad", could be directly read into BPCells in R





