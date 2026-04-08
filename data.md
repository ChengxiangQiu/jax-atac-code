---
title: Data & Code for Evolutionary Transfer Learning
---

## Count matrices

```{list-table}
- - **Cell metadata** \
    3,937,903 cells, CSV format, 616 MB
  - {button}`Download <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/df_cell.csv>`

- - **Cell x Peak matrices by** \
    h5ad format
  - {button}`Samples <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/atac_peaks_by_sample.h5ad>`
    {button}`L1 lineages <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/atac_peaks_by_L1_cell_lineage.h5ad>` \
    {button}`L2 classes <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/atac_peaks_by_L2_cell_class.h5ad>`
    {button}`L3 types <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/data_matrix/atac_peaks_by_L3_cell_type.h5ad>`

- - **Fragments** \
    by samples, L1 cell lineages, L2 cell classes, or L3 cell types
  - {button}`GEO <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE325776>`
    {button}`Script <https://github.com/ChengxiangQiu/jax-atac-code/blob/main/help_code/reorganize_fragments.py>`
```

## CREsted models

```{list-table}
- - **Evolution-naive model** \
    CREsted/1.4.0, keras/3.5.0, 62 MB
  - {button}`Download <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/CREsted_model/evolution_naive_model.tar.gz>`

- - **Evolution-aware model** \
    CREsted/1.4.0, keras/3.5.0, 67 MB
  - {button}`Download <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/CREsted_model/evolution_aware_model.tar.gz>`

- - **STEAM-v1 model** \
    CREsted/1.4.0, keras/3.13.2, 73 MB
  - {button}`Download <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/CREsted_model/STEAM_v1_model.tar.gz>`
```

## Measured & predicted tracks

```{list-table}
- - **Raw** \
    Tn5 cut sites were aggregated across cells and normalized to CPM

  - {button}`Hub <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_raw_data/hub.txt>`
    {button}`BigWig <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_raw_data/mm10/>`

- - **Evoluation-naive** \
    100-bp resolution of mm10

  - {button}`Hub <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_prediction/hub.txt>`
    {button}`BigWig <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_prediction/mm10/>`

- - **Evoluation-aware** \
    100-bp resolution of mm10

  - {button}`Hub <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_32_clusters_mouse_prediction/hub.txt>`
    {button}`BigWig <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_32_clusters_mouse_prediction/mm10/>` \
    {button}`Hub (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_32_clusters_mouse_prediction_Phred_narrow/hub.txt>`
    {button}`BigWig (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_32_clusters_mouse_prediction_Phred_narrow/mm10/>`

- - **STEAM-v1** \
    100-bp resolution of mm10

  - {button}`Hub <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_mouse_prediction/hub.txt>`
    {button}`BigWig <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_mouse_prediction/mm10/>` \
    {button}`Hub (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_mouse_prediction_Phred_narrow/hub.txt>`
    {button}`BigWig (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_mouse_prediction_Phred_narrow/mm10/>`

- - **STEAM-v1** \
    100-bp resolution of hg38

  - {button}`Hub <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_human_prediction/hub.txt>`
    {button}`BigWig <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_human_prediction/hg38/>` \
    {button}`Hub (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_human_prediction_Phred_narrow/hub.txt>`
    {button}`BigWig (GPS) <https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_human_prediction_Phred_narrow/hg38/>`
```

## Code

```{list-table}
- - **GitHub** \
    Code for data processing, analysis, and figure reproduction
  - {button}`jax-atac-code <https://github.com/ChengxiangQiu/jax-atac-code>`
```

## Contact

For questions or feedback, reach out at [Chengxiang.Qiu@dartmouth.edu](mailto:Chengxiang.Qiu@dartmouth.edu).
