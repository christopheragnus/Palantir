import palantir
import scanpy as sc
import pandas as pd
import os

import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)




ad = sc.read("./data/marrow_sample_scseq_counts.h5ad")
print(ad)

sc.pp.normalize_per_cell(ad)

palantir.preprocess.log_transform(ad)

sc.pp.highly_variable_genes(ad, n_top_genes=1500, flavor="cell_ranger")

sc.pp.pca(ad)
print(ad)

dm_res = palantir.utils.run_diffusion_maps(ad, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(ad)



terminal_states = pd.Series(
    ["DC", "Mono", "Ery"],
    index=["Run5_131097901611291", "Run5_134936662236454", "Run4_200562869397916"],
)



print(ad)


start_cell = "Run5_164698952452459"
pr_res = palantir.core.run_palantir(
    ad, start_cell, num_waypoints=500, terminal_states=terminal_states
)

