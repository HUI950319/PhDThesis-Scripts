# cisTarget Reference Databases

The following reference files are required by the SCENIC pipeline but are **not**
included in this repository because each exceeds GitHub's 100 MB file-size limit.
Download them from the Aerts Lab cisTarget resources portal and place them in
this directory before running `main.py` / `run_scenic.R`.

## Required files

### Human (hg38)

- `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (~299 MB)
- `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (~297 MB)
- `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` (~94 MB)

### Mouse (mm10)

- `mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (~227 MB)
- `mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (~226 MB)
- `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` (~108 MB)

## Download source

Official Aerts Lab cisTarget resources:

- Rankings (feather): https://resources.aertslab.org/cistarget/databases/
- Motif annotations (tbl): https://resources.aertslab.org/cistarget/motif2tf/

## Already tracked in this folder

The following small TF-list files are committed to the repo and can be used directly:

- `hsa_hgnc_tfs.motifs-v10.txt` — human TF list
- `mmu_mgi_tfs.motifs-v10.txt` — mouse TF list
