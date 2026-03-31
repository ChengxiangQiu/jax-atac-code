### Performing scGlue to integrate RNA and ATAC datasets

import anndata
import networkx as nx
import scanpy as sc
import pandas as pd
import scglue
import sys, os
from matplotlib import rcParams

work_path = ""

#############################################
### Processing two datasets, respectively ###
#############################################

### read rna dataset
rna = sc.read(os.path.join(work_path, "rna_count.mtx"))
rna.obs = pd.read_csv(os.path.join(work_path, "rna_obs.csv"), index_col = 0)
rna.var = pd.read_csv(os.path.join(work_path, "rna_var.csv"), index_col = 0)
rna.var = rna.var.drop(columns=["chr", "start", "end", "strand", "gene_ID", "gene_type", "gene_short_name"])

gene_sums = rna.X.sum(axis=0)
nonzero_genes = gene_sums > 0
rna = rna[:, nonzero_genes].copy()

rna.write(os.path.join(work_path, "rna.h5ad"), compression="gzip")

### processing rna dataset
rna.layers["counts"] = rna.X.copy()

sc.pp.highly_variable_genes(rna, n_top_genes=2500, flavor="cell_ranger")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
#sc.pp.neighbors(rna, metric="cosine")
#sc.tl.umap(rna)

### save after processing
rna.write(os.path.join(work_path, "rna_preprocessed.h5ad"), compression="gzip")

### read atac dataset
atac = sc.read(os.path.join(work_path, "atac_count.mtx"))
atac.obs = pd.read_csv(os.path.join(work_path, "atac_obs.csv"), index_col = 0)
atac.var = pd.read_csv(os.path.join(work_path, "atac_var.csv"), index_col = 0)
atac.var = atac.var.drop(columns = ["peak_ID"])

atac.write(os.path.join(work_path, "atac.h5ad"), compression="gzip")

### processing atac dataset
scglue.data.lsi(atac, n_components=100, n_iter=15)
#sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
#sc.tl.umap(atac)

### save after processing
atac.write(os.path.join(work_path, "atac_preprocessed.h5ad"), compression="gzip")


##########################
### Creating the graph ###
##########################

### obatin genomic coordinates
scglue.data.get_gene_annotation(
    rna, 
    gtf="gencode.vM12.chr_patch_hapl_scaff.annotation.update.gtf.gz", ### this can be downloaded from GENCODE (https://www.gencodegenes.org/)
    gtf_by="gene_id"
)
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

split = atac.var_names.str.split(r"[:_]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1])
atac.var["chromEnd"] = split.map(lambda x: x[2])

graph = scglue.genomics.rna_anchored_prior_graph(rna, atac)

rna.write(os.path.join(work_path, "rna_preprocessed.h5ad"), compression="gzip")
atac.write(os.path.join(work_path, "atac_preprocessed.h5ad"), compression="gzip")
nx.write_graphml(graph, os.path.join(work_path, "prior.graphml.gz"))


###################################################
### Performing scGLUE to integrate two datasets ###
###################################################

import anndata
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
import sys, os
from matplotlib import rcParams

work_Path = ""

rna = anndata.read_h5ad(os.path.join(work_path, "rna_preprocessed.h5ad"))
atac = anndata.read_h5ad(os.path.join(work_path, "atac_preprocessed.h5ad"))
graph = nx.read_graphml(os.path.join(work_path, "prior.graphml.gz"))

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))

import torch
print(torch.cuda)
print(torch.cuda.is_available())
print(torch.version.cuda)
print(torch.__version__)

glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "glue"}
)

glue.save(os.path.join(work_path, "glue.dill"))

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

import numpy as np
np.savetxt(os.path.join(work_path, "atac_X_glue.csv"), atac.obsm["X_glue"], delimiter=",")
np.savetxt(os.path.join(work_path, "rna_X_glue.csv"), rna.obsm["X_glue"], delimiter=",")

rna.obs.to_csv(os.path.join(work_path, 'rna_obs.csv'))




















