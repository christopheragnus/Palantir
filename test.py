import palantir as palantir
import dandelion as ddl
import scanpy as sc
import pandas as pd
import numpy as np
import milopy.core as milo

preprocessed_GEXdata = sc.read_h5ad("./src/data/BCR_contigs.h5ad")

# pseudobulk
GEX_adata = ddl.tl.setup_vdj_pseudobulk(preprocessed_GEXdata, mode="B",
#                                       productive_cols = ["productive_VDJ", "productive_VJ"]
                                        productive_vdj = False,
                                        productive_vj = False
                                       )                   

GEX_adata = GEX_adata[~GEX_adata.obs["celltype_annotation"].isin(["CYCLING_B", "B1", "MATURE_B"])] 


sc.pp.neighbors(GEX_adata, use_rep="X_scvi", n_neighbors=50)
milo.make_nhoods(GEX_adata, prop=0.1)

pb_GEX_adata = ddl.tl.vdj_pseudobulk(
    GEX_adata, pbs=GEX_adata.obsm["nhoods"], obs_to_take="celltype_annotation", mode="B", extract_cols=None
)

sc.tl.pca(pb_GEX_adata) # compute PCA coordinates, loadings and variance
# sc.pl.pca(pb_GEX_adata, color="celltype_annotation")
# sc.pl.pca_overview(pb_GEX_adata)

# sc.pp.neighbors(pb_GEX_adata, n_neighbors=100)
# sc.tl.umap(pb_GEX_adata)
# sc.pl.umap(pb_GEX_adata, color=["celltype_annotation"])


rootcell = np.where(pb_GEX_adata.obs["celltype_annotation"]=="LARGE_PRE_B")[0][0]
print("rootcell:", np.where(pb_GEX_adata.obs["celltype_annotation"]=="LARGE_PRE_B")[0])

terminal_states = pd.Series(
    ["IMMATURE_B"],
    index=pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="IMMATURE_B"].obs_names[
        [
            np.argmin(pb_GEX_adata[pb_GEX_adata.obs["celltype_annotation"]=="IMMATURE_B"].obsm["X_pca"][:, 1]),
        ]
    ],
)

print("terminal_states:", terminal_states)

# Run diffusion maps
pca_projections = pd.DataFrame(pb_GEX_adata.obsm["X_pca"], index=pb_GEX_adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

print(pb_GEX_adata.obs_names)

pr_res = palantir.core.run_palantir(
    ms_data,
    pb_GEX_adata.obs_names[rootcell],
    num_waypoints=500,
    terminal_states=terminal_states.index,
)

pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]
