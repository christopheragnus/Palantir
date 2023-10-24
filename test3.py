import palantir as palantir
import dandelion as ddl
import scanpy as sc
import pandas as pd
import numpy as np
import milopy.core as milo


preprocessed_GEXdata = sc.read_h5ad("./src/data/BCR_contigs.h5ad")

# subset to cells with a pair of chains
GEX_adata = ddl.tl.setup_vdj_pseudobulk(preprocessed_GEXdata, mode="B",
                                        productive_vdj = False,
                                        productive_vj = False
                                       )

GEX_adata = GEX_adata[~GEX_adata.obs["celltype_annotation"].isin(["CYCLING_B", "B1", "MATURE_B", "PLASMA_B", "LATE_PRO_B", "PRO_B", "PRE_PRO_B"])] 

sc.pp.neighbors(GEX_adata, use_rep="X_scvi", n_neighbors=50)
milo.make_nhoods(GEX_adata, prop=0.1)
sc.tl.umap(GEX_adata)

pb_GEX_adata = ddl.tl.vdj_pseudobulk(
    GEX_adata, pbs=GEX_adata.obsm["nhoods"], obs_to_take="celltype_annotation", mode="B", 
    extract_cols = None
)

# compute PCA
sc.pp.neighbors(pb_GEX_adata)
sc.tl.umap(pb_GEX_adata)

sc.tl.pca(pb_GEX_adata) # compute PCA coordinates, loadings and variance
sc.pl.pca(pb_GEX_adata, color="celltype_annotation")
# sc.pl.pca_overview(pb_GEX_adata)



# pick rootcell & terminal state
rootcell = pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="LARGE_PRE_B"].obs_names[
            np.argmin(pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="LARGE_PRE_B"].obsm["X_umap"][:, 1])
        ]

print(rootcell)

terminal_states = pd.Series(
    ["IMMATURE_B"],
    index=pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="IMMATURE_B"].obs_names[
        [
            np.argmin(pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="IMMATURE_B"].obsm["X_pca"][:, 1]),
        ]
    ],
)

print(terminal_states)

# sc.pp.neighbors(pb_GEX_adata, n_neighbors=100)

# Run diffusion maps
pca_projections = pd.DataFrame(pb_GEX_adata.obsm["X_pca"], index=pb_GEX_adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
ms_data.index = ms_data.index.astype(str)

pr_res = palantir.core.run_palantir(
    ms_data,
    rootcell,
    num_waypoints=500,
    terminal_states=terminal_states.index,
)

pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]

pb_GEX_adata = ddl.tl.pseudotime_transfer(pb_GEX_adata, pr_res)

print('plotting...')

sc.pl.pca(
    pb_GEX_adata,
    color=["pseudotime", "prob_IMMATURE_B", "celltype_annotation"],
    color_map="coolwarm"
)

bdata = ddl.tl.project_pseudotime_to_cell(
    GEX_adata, pb_GEX_adata, terminal_states.values
)
sc.pl.umap(bdata, color=["celltype_annotation"])
sc.pl.umap(
    bdata,
    color=["pseudotime", "prob_IMMATURE_B"],
    color_map="coolwarm",
)